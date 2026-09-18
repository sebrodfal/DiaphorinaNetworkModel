(* ::Package:: *)

(* Network experiment infrastructure for DiaphorinaModel.
   Reuses the ODE structure defined in SystemEvolveNetwork (Functions.m) but
   decouples solving from plotting, so batches of experiments can run headlessly
   (via wolframscript, no .nb generated) and results can be exported as data. *)


(* ---------- Solver (Model + Solver, no plotting) ---------- *)

ClearAll[SolveNetworkData];
SolveNetworkData[graph_?GraphQ, Dd_, Dt_, \[Delta]_, T_,
    {{r1_, b12_, b13_, a1_, c1_},
     {r2_, b21_, b24_, a2_, c2_},
     {r3_, b31_, a3_, c3_},
     {r4_, b42_, a4_, c4_}},
    {Di0_, B0_, x30_, x40_, tt_},
    releaseNodes_: All
] := Module[{L, n, eqs, ics, vars, sol, events, targets},
    L = Normal @ KirchhoffMatrix[graph];
    n = VertexCount[graph];

    eqs = Flatten @ {
        Table[
            x1[i]'[t] == x1[i][t]*(r1 + b12*x2[i][t] + b13*x3[i][t] -
                (a1 + c1*(b12*x2[i][t] + b13*x3[i][t]))*x1[i][t]) -
                Dd*Sum[L[[i, j]] x1[j][t], {j, n}],
            {i, n}
        ],
        Table[
            x2[i]'[t] == x2[i][t]*(r2 + b21*x1[i][t] + b24*x4[i][t] -
                (a2 + c2*(b21*x1[i][t] + b24*x4[i][t]))*x2[i][t]),
            {i, n}
        ],
        Table[
            x3[i]'[t] == x3[i][t]*(r3 + b31*x1[i][t] -
                (a3 + c3*b31*x1[i][t])*x3[i][t]) -
                Dt*Sum[L[[i, j]] x3[j][t], {j, n}],
            {i, n}
        ],
        Table[
            x4[i]'[t] == x4[i][t]*(r4 + b42*x2[i][t] -
                (a4 + c4*b42*x2[i][t])*x4[i][t]),
            {i, n}
        ]
    };

    ics = Flatten @ {
        Table[x1[i][0] == Di0[[i]], {i, n}],
        Table[x2[i][0] == B0[[i]], {i, n}],
        Table[x3[i][0] == x30[[i]], {i, n}],
        Table[x4[i][0] == x40[[i]], {i, n}]
    };

    vars = Flatten @ {
        Table[x1[i], {i, n}],
        Table[x2[i], {i, n}],
        Table[x3[i], {i, n}],
        Table[x4[i], {i, n}]
    };

    targets = If[releaseNodes === All, Range[n], releaseNodes];
    events = If[\[Delta] == 0 || targets === {}, {},
        {WhenEvent[Mod[t, T] == 0 && t > 0,
            Evaluate @ Table[x3[i][t] -> x3[i][t] + \[Delta], {i, targets}]]}
    ];

    sol = First @ NDSolve[
        Join[eqs, ics, events],
        vars, {t, 0, tt},
        Method -> {"TimeIntegration" -> "ExplicitRungeKutta"}
    ];

    <|
        "Solution" -> sol, "Graph" -> graph, "N" -> n,
        "Dd" -> Dd, "Dt" -> Dt, "Delta" -> \[Delta], "T" -> T, "tt" -> tt,
        "ReleaseNodes" -> targets
    |>
];


(* ---------- Topology generators ---------- *)

ClearAll[MakeTopology];
MakeTopology[name_String, n_Integer, seed_: 1] := (
    SeedRandom[seed];
    Switch[name,
        "Path", PathGraph[Range[n]],
        "Cycle", CycleGraph[n],
        "Star", StarGraph[n],
        "Complete", CompleteGraph[n],
        "Random", RandomGraph[BernoulliGraphDistribution[n, 0.4]],
        "SmallWorld", RandomGraph[WattsStrogatzGraphDistribution[n, 0.2, Max[1, Min[2, Floor[n/2] - 1]]]],
        "ScaleFree", RandomGraph[BarabasiAlbertGraphDistribution[n, Min[2, n - 1]]],
        _, Missing["UnknownTopology", name]
    ]
);


(* ---------- Metrics ---------- *)

ClearAll[SampleTimes];
SampleTimes[data_Association, npts_: 400] := Subdivide[0., data["tt"], npts];

ClearAll[NodeIF];
NodeIF[data_Association, var_, i_] := var[i] /. data["Solution"];

ClearAll[AUC];
AUC[data_Association, var_, i_] := NIntegrate[NodeIF[data, var, i][t], {t, 0, data["tt"]}];

ClearAll[TotalAUC];
TotalAUC[data_Association, var_] := Sum[AUC[data, var, i], {i, data["N"]}];

ClearAll[LastCycleMean];
LastCycleMean[data_Association, var_, i_] := Module[{T = data["T"], tt = data["tt"], t0, if},
    if = NodeIF[data, var, i];
    If[T <= 0 || T >= tt, Return[if[tt]]];
    t0 = tt - Mod[tt, T] - T;
    If[t0 < 0, t0 = 0];
    NIntegrate[if[t], {t, t0, tt}]/(tt - t0)
];

ClearAll[WindowMean];
(* Mean over an explicit [t0,t1] window, independent of the release period T. Exists
   alongside LastCycleMean (which averages over the single most recent release cycle)
   so results here can be compared directly against the single-node paper's own
   reporting convention -- Table 1 / ReleaseSweep / ScenarioMetrics in ReproduceResults.m
   all report the mean over t in [tt/2,tt], not over one release cycle, and with T=15
   on a tt=100 horizon those are very different windows (~15 units vs ~50 units). *)
WindowMean[data_Association, var_, i_, win_List] := Module[{if},
    if = NodeIF[data, var, i];
    NIntegrate[if[t], {t, win[[1]], win[[2]]}]/(win[[2]] - win[[1]])
];

ClearAll[PeakValue];
PeakValue[data_Association, var_, i_] := Max[NodeIF[data, var, i] /@ SampleTimes[data]];

ClearAll[NodeMinValue];
NodeMinValue[data_Association, var_, i_] := Min[NodeIF[data, var, i] /@ SampleTimes[data]];

ClearAll[SpatialCV];
SpatialCV[data_Association, var_] := Module[{vals},
    vals = Table[LastCycleMean[data, var, i], {i, data["N"]}];
    If[Mean[vals] == 0, 0., StandardDeviation[vals]/Mean[vals]]
];

ClearAll[SyncIndex];
SyncIndex[data_Association, var_] := Module[{n = data["N"], tvals, series, pairs, corrs},
    If[n < 2, Return[1.]];
    tvals = SampleTimes[data];
    series = Table[NodeIF[data, var, i] /@ tvals, {i, n}];
    pairs = Subsets[Range[n], {2}];
    corrs = Correlation[series[[#[[1]]]], series[[#[[2]]]]] & /@ pairs;
    N[Mean[corrs]]
];

ClearAll[AlgebraicConnectivity];
AlgebraicConnectivity[graph_?GraphQ] := Module[{ev},
    ev = Sort[Eigenvalues[N @ Normal @ KirchhoffMatrix[graph]]];
    If[Length[ev] < 2, 0., ev[[2]]]
];

ClearAll[ExperimentMetrics];
Options[ExperimentMetrics] = {"Window" -> Automatic};
(* "Window" defaults to {tt/2,tt}, matching ReleaseSweep/ScenarioMetrics in the
   single-node ReproduceResults.m, so "PestWindowMean"/"ParasitoidWindowMean" here are
   directly comparable to that paper's MeanDensity numbers (e.g. the network's own
   disconnected-topology limit should reproduce Table 1's 91.3/34.2/21.9). *)
ExperimentMetrics[data_Association, opts:OptionsPattern[]] := Module[{win},
    win = OptionValue["Window"];
    If[win === Automatic, win = {data["tt"]/2, data["tt"]}];
    <|
        "PestPeak" -> Max[Table[PeakValue[data, x1, i], {i, data["N"]}]],
        "PestMin" -> Min[Table[NodeMinValue[data, x1, i], {i, data["N"]}]],
        "PestAUC" -> TotalAUC[data, x1],
        "PestLastCycleMean" -> Mean[Table[LastCycleMean[data, x1, i], {i, data["N"]}]],
        "PestWindowMean" -> Mean[Table[WindowMean[data, x1, i, win], {i, data["N"]}]],
        "PestSpatialCV" -> SpatialCV[data, x1],
        "PestSync" -> SyncIndex[data, x1],
        "ParasitoidLastCycleMean" -> Mean[Table[LastCycleMean[data, x3, i], {i, data["N"]}]],
        "ParasitoidWindowMean" -> Mean[Table[WindowMean[data, x3, i, win], {i, data["N"]}]],
        "ParasitoidSync" -> SyncIndex[data, x3],
        "AlgebraicConnectivity" -> AlgebraicConnectivity[data["Graph"]],
        "MeanDegree" -> N @ Mean[VertexDegree[data["Graph"]]]
    |>
];
