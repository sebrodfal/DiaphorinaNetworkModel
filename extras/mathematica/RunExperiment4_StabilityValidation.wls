#!/usr/bin/env wolframscript
(* Experiment 4: validates the new Master-Stability-without-periodicity machinery
   (NetworkStability.wl's MasterLyapunovExponent/MasterLyapunovSpectrum) and applies
   it to the real topologies used in Experiments 2/3.

   Part A -- ground truth: mu=0 must reduce to the plain single-node problem,
   independent of Dd/Dt (mu*Dmat vanishes identically at mu=0). Validated two ways:
   (i) against an independently-written raw single-node Lyapunov calculation that
   does not go through ReferenceTrajectory/PerturbationSegment at all, and (ii)
   internally, by checking MasterLyapunovExponent gives the same value at mu=0
   regardless of which Dd/Dt is passed in.

   Part B -- applies MasterLyapunovSpectrum to the ACTUAL Laplacian spectrum of
   each of the 7 topologies from Experiments 2/3 (n=8, baseline delta=35,T=15,
   Dd=Dt=0.03), predicting whether the synchronized regime is expected to be
   stable (every transverse mode's Lyapunov exponent < 0). Cross-references that
   prediction against the empirical PestSync/PestSpatialCV already measured for
   the same (Dd,Dt,Strategy=Uniform) slice in Experiment 3's summary CSV -- this
   is the first place theory (graph spectrum) and simulation actually meet.

   Runs entirely headless, no .nb generated. *)

SetDirectory[DirectoryName[$InputFileName]];
Get["NetworkExperiments.wl"]; (* for MakeTopology *)
Get["NetworkStability.wl"];
If[! DirectoryQ["Results"], CreateDirectory["Results"]];

p1 = {-0.20, 0.040, -0.010, 5.55*10^-4, 0.005};
p2 = {-0.72, -0.010, 0.054, 5.55*10^-4, 0.005};
p3 = {-0.30, 0.010, 5.55*10^-4, 0.005};
p4 = {0.72, -0.072, 5.55*10^-4, 0.005};
params = {p1, p2, p3, p4};

deltaBase = 35; Tbase = 15; DdBase = 0.03; DtBase = 0.03; n = 8;


(* ============================================================================
   PART A -- ground truth check: mu=0 must reduce to the single-node problem
   ============================================================================ *)
Print["### PART A: mu=0 reduction check ###"];

(* Independent (from-scratch) single-node Lyapunov exponent -- does not call
   ReferenceTrajectory or PerturbationSegment, so it exercises none of the same
   code as MasterLyapunovExponent; only LocalField/LocalJacobianFunction are
   shared, making this a genuine independent cross-check. *)
ClearAll[RawSingleNodeLyapunov];
RawSingleNodeLyapunov[params_, \[Delta]_, T_, tTotal_, burnFrac_, segLen_] := Module[
    {f, Jf, events, sSol, sFns, burn, nSeg, step, states},
    f[{y1_, y2_, y3_, y4_}] := LocalField[params][{y1, y2, y3, y4}];
    Jf = LocalJacobianFunction[params];
    events = If[\[Delta] == 0, {}, {WhenEvent[Mod[t, T] == 0 && t > 0, w3[t] -> w3[t] + \[Delta]]}];
    sSol = First @ NDSolve[
        Join[{
            w1'[t] == f[{w1[t], w2[t], w3[t], w4[t]}][[1]],
            w2'[t] == f[{w1[t], w2[t], w3[t], w4[t]}][[2]],
            w3'[t] == f[{w1[t], w2[t], w3[t], w4[t]}][[3]],
            w4'[t] == f[{w1[t], w2[t], w3[t], w4[t]}][[4]],
            w1[0] == 50, w2[0] == 10, w3[0] == 15, w4[0] == 30
        }, events], {w1, w2, w3, w4}, {t, 0, tTotal},
        Method -> {"TimeIntegration" -> "ExplicitRungeKutta"}
    ];
    sFns = {w1, w2, w3, w4} /. sSol;
    burn = burnFrac*tTotal;
    nSeg = Floor[(tTotal - burn)/segLen];
    step[{tC_, vC_, logSum_}] := Module[{s1f, s2f, s3f, s4f, eqs, ics, sol, vEnd, nrm},
        {s1f, s2f, s3f, s4f} = sFns;
        eqs = Table[u[k]'[t] == Sum[Jf[s1f[t], s2f[t], s3f[t], s4f[t]][[k, j]]*u[j][t], {j, 1, 4}], {k, 1, 4}];
        ics = Table[u[k][tC] == vC[[k]], {k, 1, 4}];
        sol = First @ NDSolve[Join[eqs, ics], Table[u[k], {k, 1, 4}], {t, tC, tC + segLen},
            Method -> {"TimeIntegration" -> "ExplicitRungeKutta"}];
        vEnd = Table[(u[k] /. sol)[tC + segLen], {k, 1, 4}];
        nrm = Norm[vEnd];
        {tC + segLen, vEnd/nrm, logSum + Log[nrm]}
    ];
    states = NestList[step, {burn, {1., 0., 0., 0.}, 0.}, nSeg];
    Last[states[[All, 3]]]/(nSeg*segLen)
];

rawLambda0 = RawSingleNodeLyapunov[params, deltaBase, Tbase, 1200, 0.3, Tbase];
viaMasterDd0 = MasterLyapunovExponent[params, 0, 0, 0, deltaBase, Tbase]["LyapunovExponent"];
viaMasterDdReal = MasterLyapunovExponent[params, DdBase, DtBase, 0, deltaBase, Tbase]["LyapunovExponent"];

Print["Raw single-node lambda (independent calc):      ", rawLambda0];
Print["MasterLyapunovExponent, mu=0, Dd=Dt=0:           ", viaMasterDd0];
Print["MasterLyapunovExponent, mu=0, Dd=Dt=", DdBase, ":  ", viaMasterDdReal];
Print["Max pairwise diff (should be ~0): ",
    Max[Abs[{rawLambda0, viaMasterDd0, viaMasterDdReal} - rawLambda0]]];
partAPass = Max[Abs[{rawLambda0, viaMasterDd0, viaMasterDdReal} - rawLambda0]] < 10^-4;
Print[If[partAPass, "PASS: mu=0 reduces correctly to the single-node problem, independent of Dd/Dt.",
    "FAIL: mismatch -- check ReferenceTrajectory/PerturbationSegment against RawSingleNodeLyapunov."]];


(* ============================================================================
   PART B -- MasterLyapunovSpectrum on the real topologies from Experiments 2/3
   ============================================================================ *)
Print["\n### PART B: synchronized-state stability prediction per topology ###"];

topologies = {"Path", "Cycle", "Star", "Complete", "Random", "SmallWorld", "ScaleFree"};

topologyResult[topo_] := Module[{graph, muAll, muTransverse, spectrum, maxLambda, reliableAll},
    graph = MakeTopology[topo, n];
    muAll = Sort[Eigenvalues[N @ Normal @ KirchhoffMatrix[graph]]];
    muTransverse = Rest[muAll]; (* drop mu_1=0, the synchronized mode itself *)
    spectrum = MasterLyapunovSpectrum[params, DdBase, DtBase, muTransverse, deltaBase, Tbase];
    maxLambda = Max[spectrum[[All, "LyapunovExponent"]]];
    reliableAll = And @@ spectrum[[All, "Reliable"]];
    <|
        "Topology" -> topo, "AlgebraicConnectivity" -> muAll[[2]],
        "MaxTransverseLyapunov" -> maxLambda,
        "PredictedSyncStable" -> (maxLambda < 0),
        "AllReliable" -> reliableAll
    |>
];

results = topologyResult /@ topologies;

csvPath = FileNameJoin[{"Results", "experiment4_stability_prediction.csv"}];
columnNames = Keys[results[[1]]];
csvRows = Prepend[Values /@ results, columnNames];
Export[csvPath, csvRows, "CSV"];
Print["Saved: ", csvPath];

Scan[
    Print["  ", #Topology, ": max transverse lambda=", #MaxTransverseLyapunov,
        "  predicted stable sync=", #PredictedSyncStable, "  all segments reliable=", #AllReliable] &,
    results
];

(* Cross-reference against the empirical PestSync/PestSpatialCV already measured in
   Experiment 3 (same Dd=Dt=0.03, Uniform strategy) -- reads the prior run's summary
   CSV directly since Experiment 3 has not been re-run inside this project yet. *)
oldSummaryPath = "C:\\Users\\srodrifa\\Documents\\Model (depricated)\\Model (depricated)\\Test\\Results\\experiment3_replicates_summary.csv";
If[FileExistsQ[oldSummaryPath],
    Module[{summary, bySlice, comparison},
        summary = Import[oldSummaryPath, "Dataset", HeaderLines -> 1];
        bySlice = Normal[summary[Select[Abs[#Dd - DdBase] < 10^-9 && Abs[#Dt - DtBase] < 10^-9 && #Strategy == "Uniform" &]]];
        comparison = Table[
            Module[{pred = SelectFirst[results, #Topology == topo &], emp = SelectFirst[bySlice, #Topology == topo &]},
                <|"Topology" -> topo, "PredictedSyncStable" -> pred["PredictedSyncStable"],
                  "MaxTransverseLyapunov" -> pred["MaxTransverseLyapunov"],
                  "EmpiricalPestSync_Mean" -> emp["PestSync_Mean"],
                  "EmpiricalPestSpatialCV_Mean" -> emp["PestSpatialCV_Mean"]|>
            ],
            {topo, topologies}
        ];
        Print["\n### Theory vs. simulation (Experiment 3, Dd=Dt=0.03, Uniform) ###"];
        Scan[
            Print["  ", #Topology, ": predicted stable=", #PredictedSyncStable,
                "  (lambda=", #MaxTransverseLyapunov, ")   empirical PestSync=", #["EmpiricalPestSync_Mean"],
                "  PestSpatialCV=", #["EmpiricalPestSpatialCV_Mean"]] &,
            comparison
        ];
        compPath = FileNameJoin[{"Results", "experiment4_theory_vs_simulation.csv"}];
        compColumnNames = Keys[comparison[[1]]];
        compCsvRows = Prepend[Values /@ comparison, compColumnNames];
        Export[compPath, compCsvRows, "CSV"];
        Print["Saved: ", compPath];
    ],
    Print["\nSkipped theory-vs-simulation cross-reference: could not find ", oldSummaryPath];
];

Print["\nDONE."];
