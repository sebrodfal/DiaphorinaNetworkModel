#!/usr/bin/env wolframscript
(* Experiment 3: replicate stochastic topologies to quantify variance.
   "Random", "SmallWorld", "ScaleFree" are each a single graph *instance* drawn
   from a distribution -- Experiment 2 used one fixed seed per topology, so any
   difference attributed to "Random" vs "Path" was confounded with which
   particular random graph got drawn. Here we redraw each stochastic topology
   nSeeds times (independent graph instances) and re-run the same point-source
   outbreak scenario as Experiment 2, so downstream analysis can report
   mean +/- std per (topology, Dd, Dt, strategy) instead of a single sample.

   Deterministic topologies (Path, Cycle, Star, Complete) are not redrawn --
   they are the same graph regardless of seed -- but are included once for
   reference/comparison in the same output files.

   Runs entirely headless (wolframscript -file), no .nb is created. *)

SetDirectory[DirectoryName[$InputFileName]];
Get["NetworkExperiments.wl"];
If[! DirectoryQ["Results"], CreateDirectory["Results"]];

p1 = {-0.20, 0.040, -0.010, 5.55*10^-4, 0.005};
p2 = {-0.72, -0.010, 0.054, 5.55*10^-4, 0.005};
p3 = {-0.30, 0.010, 5.55*10^-4, 0.005};
p4 = {0.72, -0.072, 5.55*10^-4, 0.005};
params = {p1, p2, p3, p4};

n = 8;
tt = 100;
T = 15;
deltaAmt = 35;

ics = {
    Prepend[ConstantArray[20, n - 1], 200],
    ConstantArray[10, n],
    ConstantArray[15, n],
    ConstantArray[30, n],
    tt
};

deterministicTopologies = {"Path", "Cycle", "Star", "Complete"};
stochasticTopologies = {"Random", "SmallWorld", "ScaleFree"};
nSeeds = 10;
diffusionGrid = {0, 0.01, 0.03, 0.1};
strategies = {"None", "Uniform", "TargetedAtSource"};

results = {};
runCount = 0;

runOne[topo_, seed_, graph_] := Module[{},
    Do[
        Do[
            {\[Delta], releaseNodes, strategyLabel} = Switch[strategy,
                "None", {0, {}, "None"},
                "Uniform", {deltaAmt, All, "Uniform"},
                "TargetedAtSource", {deltaAmt, {1}, "TargetedAtSource"}
            ];
            data = SolveNetworkData[graph, Dd, Dt, \[Delta], T, params, ics, releaseNodes];
            metrics = ExperimentMetrics[data];
            AppendTo[results, Join[
                <|
                    "Topology" -> topo, "Seed" -> seed, "N" -> n, "Dd" -> Dd, "Dt" -> Dt,
                    "Strategy" -> strategyLabel, "Delta" -> \[Delta], "T" -> T
                |>,
                metrics
            ]];
            runCount++;
            ,
            {strategy, strategies}
        ],
        {Dd, diffusionGrid}, {Dt, diffusionGrid}
    ];
];

(* Deterministic topologies: single instance, Seed -> 1 (graph does not depend on seed) *)
Do[
    graph = MakeTopology[topo, n, 1];
    runOne[topo, 1, graph],
    {topo, deterministicTopologies}
];

(* Stochastic topologies: nSeeds independent graph draws each *)
Do[
    Do[
        graph = MakeTopology[topo, n, seed];
        runOne[topo, seed, graph],
        {seed, nSeeds}
    ],
    {topo, stochasticTopologies}
];

Print["Completed ", runCount, " runs."];

(* Raw, replicate-level CSV *)
rawPath = FileNameJoin[{"Results", "experiment3_replicates_raw.csv"}];
columnNames = Keys[results[[1]]];
csvRows = Prepend[Values /@ results, columnNames];
Export[rawPath, csvRows, "CSV"];
Print["Saved: ", rawPath];

(* Summary CSV: mean/std per (Topology, Dd, Dt, Strategy), collapsing over Seed *)
ds = Dataset[results];
metricKeys = {"PestPeak", "PestMin", "PestAUC", "PestLastCycleMean", "PestSpatialCV",
    "PestSync", "ParasitoidLastCycleMean", "ParasitoidSync", "AlgebraicConnectivity", "MeanDegree"};

grouped = GroupBy[results, {#Topology, #Dd, #Dt, #Strategy} &];

summaryRows = KeyValueMap[
    Function[{key, group},
        Join[
            <|"Topology" -> key[[1]], "Dd" -> key[[2]], "Dt" -> key[[3]], "Strategy" -> key[[4]],
              "NumSeeds" -> Length[group]|>,
            Association @ Flatten @ Table[
                Module[{vals = Lookup[group, m]},
                    {
                        (m <> "_Mean") -> N[Mean[vals]],
                        (m <> "_Std") -> If[Length[vals] > 1, N[StandardDeviation[vals]], 0.]
                    }
                ],
                {m, metricKeys}
            ]
        ]
    ],
    grouped
];

summaryPath = FileNameJoin[{"Results", "experiment3_replicates_summary.csv"}];
summaryColumnNames = Keys[summaryRows[[1]]];
summaryCsvRows = Prepend[Values /@ summaryRows, summaryColumnNames];
Export[summaryPath, summaryCsvRows, "CSV"];
Print["Saved: ", summaryPath];

(* Figure: PestAUC mean +/- std by topology, Strategy="Uniform", Dd=Dt=0.03 *)
allTopologies = Join[deterministicTopologies, stochasticTopologies];
slice = Select[summaryRows, #Dd == 0.03 && #Dt == 0.03 && #Strategy == "Uniform" &];
sliceByTopo = Association[#Topology -> # & /@ slice];

meanVals = Table[sliceByTopo[topo]["PestAUC_Mean"], {topo, allTopologies}];
stdVals = Table[sliceByTopo[topo]["PestAUC_Std"], {topo, allTopologies}];

errorPlot = ErrorListPlot[
    Transpose[{Range[Length[allTopologies]], meanVals, stdVals}] /. {i_, m_, s_} -> {i, m} \[PlusMinus] s,
    PlotStyle -> PointSize[0.02],
    Frame -> True,
    FrameLabel -> {"Topology", "Pest AUC (network-wide, Uniform release, Dd=Dt=0.03)"},
    FrameTicks -> {{Automatic, None}, {Table[{i, allTopologies[[i]]}, {i, Length[allTopologies]}], None}},
    GridLines -> Automatic, GridLinesStyle -> Directive[Dashed],
    PlotLabel -> "Mean \[PlusMinus] std of Pest AUC across " <> ToString[nSeeds] <> " seeds (stochastic topologies)",
    ImageSize -> Large
];
pngPath = FileNameJoin[{"Results", "experiment3_summary_errorbars.png"}];
Export[pngPath, errorPlot];
Print["Saved: ", pngPath];
