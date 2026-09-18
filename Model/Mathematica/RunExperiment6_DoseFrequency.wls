#!/usr/bin/env wolframscript
(* Experiment 6: dose-frequency sweep (delta x T) crossed with topology, over
   the SAME short/medium horizon used throughout the preprint (tt=100). We
   deliberately stay inside this horizon: Experiment 5's bifurcation scan
   showed that both the single-node model AND the raw continuous system
   (delta=0) develop wide, irregular oscillations over long horizons
   (hundreds of release cycles) -- consistent only over ~100 time units,
   matching what the preprint and notebooks actually demonstrate. So this
   sweep asks a horizon-appropriate question: within the timeframe the model
   is validated for, which (delta,T) combination controls the point-source
   outbreak most cost-effectively, and does that answer depend on topology?

   Cost-effectiveness: TotalReleased = delta * (number of release events in
   [0,tt]) * (number of release nodes). CostEfficiency = AUC avoided per unit
   of parasitoid released, relative to the no-control baseline for the same
   topology/diffusion.

   Runs entirely headless, no .nb generated. *)

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
DdVal = 0.03;
DtVal = 0.03;

ics = {
    Prepend[ConstantArray[20, n - 1], 200],
    ConstantArray[10, n],
    ConstantArray[15, n],
    ConstantArray[30, n],
    tt
};

topologies = {"Path", "Star", "Complete", "Random"};
deltaGrid = {10, 20, 30, 40, 50};
Tgrid = {5, 10, 15, 20, 30};
strategies = {"Uniform", "TargetedAtSource"};

results = {};
baselineAUC = <||>;

Do[
    graph = MakeTopology[topo, n];

    (* No-control baseline for this topology, to normalize cost-effectiveness *)
    baseData = SolveNetworkData[graph, DdVal, DtVal, 0, 15, params, ics, {}];
    baselineAUC[topo] = ExperimentMetrics[baseData]["PestAUC"];

    Do[
        Do[
            Do[
                releaseNodes = If[strategy == "Uniform", All, {1}];
                nReleaseNodes = If[strategy == "Uniform", n, 1];
                data = SolveNetworkData[graph, DdVal, DtVal, deltaAmt, Tval, params, ics, releaseNodes];
                metrics = ExperimentMetrics[data];
                nEvents = Floor[tt/Tval];
                totalReleased = deltaAmt*nEvents*nReleaseNodes;
                aucAvoided = baselineAUC[topo] - metrics["PestAUC"];
                AppendTo[results, Join[
                    <|
                        "Topology" -> topo, "Delta" -> deltaAmt, "T" -> Tval, "Strategy" -> strategy,
                        "TotalReleased" -> totalReleased, "BaselineAUC" -> baselineAUC[topo],
                        "AUCAvoided" -> aucAvoided,
                        "CostEfficiency" -> If[totalReleased > 0, aucAvoided/totalReleased, 0]
                    |>,
                    metrics
                ]];
                ,
                {strategy, strategies}
            ],
            {Tval, Tgrid}
        ],
        {deltaAmt, deltaGrid}
    ],
    {topo, topologies}
];

Print["Completed ", Length[results], " runs."];

csvPath = FileNameJoin[{"Results", "experiment6_dose_frequency.csv"}];
columnNames = Keys[results[[1]]];
csvRows = Prepend[Values /@ results, columnNames];
Export[csvPath, csvRows, "CSV"];
Print["Saved: ", csvPath];

(* Best (delta,T) per topology x strategy by CostEfficiency *)
grouped = GroupBy[results, {#Topology, #Strategy} &];
best = KeyValueMap[
    Function[{key, group}, MaximalBy[group, #CostEfficiency &][[1]]],
    grouped
];
bestPath = FileNameJoin[{"Results", "experiment6_best_dose_frequency.csv"}];
bestColumnNames = Keys[best[[1]]];
bestCsvRows = Prepend[Values /@ best, bestColumnNames];
Export[bestPath, bestCsvRows, "CSV"];
Print["Saved: ", bestPath];
Do[
    Print[b["Topology"], " / ", b["Strategy"], ": best delta=", b["Delta"], " T=", b["T"],
        "  CostEfficiency=", b["CostEfficiency"], "  PestAUC=", b["PestAUC"]],
    {b, best}
];
