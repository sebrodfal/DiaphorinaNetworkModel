#!/usr/bin/env wolframscript
(* Experiment 2: point-source outbreak containment across network topologies.
   Node 1 starts as an outbreak focus (high Diaphorina); all other nodes start
   pest-free-ish. We sweep topology, pest/parasitoid diffusion rates (Dd,Dt),
   and release strategy (no control / uniform release at every node / targeted
   release only at the outbreak node) to see whether network structure and
   parasitoid dispersal can substitute for network-wide releases.

   Runs entirely headless (wolframscript -file), no .nb is created.
   Output: a long-format CSV in Results/, one row per run. *)

SetDirectory[DirectoryName[$InputFileName]];
Get["NetworkExperiments.wl"];
If[! DirectoryQ["Results"], CreateDirectory["Results"]];

(* Base parameters, from tabla_latex.txt (single-node calibration in the preprint) *)
p1 = {-0.20, 0.040, -0.010, 5.55*10^-4, 0.005};
p2 = {-0.72, -0.010, 0.054, 5.55*10^-4, 0.005};
p3 = {-0.30, 0.010, 5.55*10^-4, 0.005};
p4 = {0.72, -0.072, 5.55*10^-4, 0.005};
params = {p1, p2, p3, p4};

n = 8;
tt = 100;
T = 15;
deltaAmt = 35;

(* Point-source outbreak: node 1 starts at 200 (outbreak), rest at 20 (baseline) *)
ics = {
    Prepend[ConstantArray[20, n - 1], 200],
    ConstantArray[10, n],
    ConstantArray[15, n],
    ConstantArray[30, n],
    tt
};

topologies = {"Path", "Cycle", "Star", "Complete", "Random", "SmallWorld", "ScaleFree"};
diffusionGrid = {0, 0.01, 0.03, 0.1};
strategies = {"None", "Uniform", "TargetedAtSource"};

results = {};
runCount = 0;

Do[
    graph = MakeTopology[topo, n];
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
                    "Topology" -> topo, "N" -> n, "Dd" -> Dd, "Dt" -> Dt,
                    "Strategy" -> strategyLabel, "Delta" -> \[Delta], "T" -> T
                |>,
                metrics
            ]];
            runCount++;
            ,
            {strategy, strategies}
        ],
        {Dd, diffusionGrid}, {Dt, diffusionGrid}
    ],
    {topo, topologies}
];

Print["Completed ", runCount, " runs."];

ds = Dataset[results];
csvPath = FileNameJoin[{"Results", "experiment2_outbreak_containment.csv"}];
columnNames = Keys[results[[1]]];
csvRows = Prepend[Values /@ results, columnNames];
Export[csvPath, csvRows, "CSV"];
Print["Saved: ", csvPath];

(* Quick summary figure: mean PestAUC by topology x strategy, at Dd=Dt=0.03
   (a representative mid-range dispersal rate), to sanity-check the pipeline
   end-to-end without opening a notebook. *)
subset = Normal[ds[Select[#Dd == 0.03 && #Dt == 0.03 &]]];
grouped = GroupBy[subset, #Strategy &, Mean[#PestAUC & /@ #] &];

barChart = BarChart[
    Table[Lookup[grouped, s, 0], {s, strategies}],
    ChartLabels -> strategies,
    ChartLegends -> None,
    PlotLabel -> "Mean network Pest AUC by release strategy (Dd=Dt=0.03, avg over topologies)",
    FrameLabel -> {"Release strategy", "Total Diaphorina AUC (network-wide)"},
    Frame -> True, GridLines -> Automatic, GridLinesStyle -> Directive[Dashed],
    ImageSize -> Large
];
pngPath = FileNameJoin[{"Results", "experiment2_summary_strategy.png"}];
Export[pngPath, barChart];
Print["Saved: ", pngPath];
