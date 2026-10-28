#!/usr/bin/env wolframscript
(* Experiment 7: one-at-a-time local sensitivity analysis of the 18 local
   parameters (r,b,a,c for the four species) on the point-source outbreak
   scenario, within the validated short horizon (tt=100). Each parameter is
   perturbed +/-10% from its tabla_latex.txt value, holding everything else
   fixed (topology=Complete, Dd=Dt=0.03, delta=35, T=15, Uniform release),
   and we record the resulting % change in PestAUC and PestLastCycleMean.
   This ranks which biological/ecological parameters the control outcome is
   most sensitive to -- useful both as a robustness check and to prioritize
   which parameters would most need field calibration.

   Runs entirely headless, no .nb generated. *)

SetDirectory[DirectoryName[$InputFileName]];
Get["NetworkExperiments.wl"];
If[! DirectoryQ["Results"], CreateDirectory["Results"]];

baseP1 = {-0.20, 0.040, -0.010, 5.55*10^-4, 0.005};
baseP2 = {-0.72, -0.010, 0.054, 5.55*10^-4, 0.005};
baseP3 = {-0.30, 0.010, 5.55*10^-4, 0.005};
baseP4 = {0.72, -0.072, 5.55*10^-4, 0.005};
baseParams = {baseP1, baseP2, baseP3, baseP4};

paramNames = {
    {1, 1, "r1"}, {1, 2, "b12"}, {1, 3, "b13"}, {1, 4, "a1"}, {1, 5, "c1"},
    {2, 1, "r2"}, {2, 2, "b21"}, {2, 3, "b24"}, {2, 4, "a2"}, {2, 5, "c2"},
    {3, 1, "r3"}, {3, 2, "b31"}, {3, 3, "a3"}, {3, 4, "c3"},
    {4, 1, "r4"}, {4, 2, "b42"}, {4, 3, "a4"}, {4, 4, "c4"}
};

n = 8;
tt = 100;
DdVal = 0.03; DtVal = 0.03;
deltaAmt = 35; Tval = 15;

ics = {
    Prepend[ConstantArray[20, n - 1], 200],
    ConstantArray[10, n], ConstantArray[15, n], ConstantArray[30, n], tt
};

graph = MakeTopology["Complete", n];

RunWithParams[p_] := ExperimentMetrics[SolveNetworkData[graph, DdVal, DtVal, deltaAmt, Tval, p, ics, All]];

baseMetrics = RunWithParams[baseParams];
Print["Baseline PestAUC = ", baseMetrics["PestAUC"], "   PestLastCycleMean = ", baseMetrics["PestLastCycleMean"]];

Perturb[block_, idx_, factor_] := ReplacePart[baseParams, {block, idx} -> baseParams[[block, idx]]*factor];

results = {};
Do[
    {block, idx, name} = pn;
    Do[
        pParams = Perturb[block, idx, factor];
        m = RunWithParams[pParams];
        AppendTo[results, <|
            "Parameter" -> name, "Factor" -> factor,
            "PerturbedValue" -> baseParams[[block, idx]]*factor,
            "PestAUC" -> m["PestAUC"],
            "PctChangeAUC" -> 100.*(m["PestAUC"] - baseMetrics["PestAUC"])/baseMetrics["PestAUC"],
            "PestLastCycleMean" -> m["PestLastCycleMean"],
            "PctChangeLastCycleMean" -> 100.*(m["PestLastCycleMean"] - baseMetrics["PestLastCycleMean"])/baseMetrics["PestLastCycleMean"]
        |>];
        ,
        {factor, {0.9, 1.1}}
    ],
    {pn, paramNames}
];

Print["Completed ", Length[results], " perturbation runs."];

csvPath = FileNameJoin[{"Results", "experiment7_sensitivity.csv"}];
columnNames = Keys[results[[1]]];
csvRows = Prepend[Values /@ results, columnNames];
Export[csvPath, csvRows, "CSV"];
Print["Saved: ", csvPath];

(* Sensitivity score per parameter: max |%change in AUC| across the two perturbations *)
grouped = GroupBy[results, #Parameter &];
sensitivity = KeyValueMap[
    <|"Parameter" -> #1, "MaxAbsPctChangeAUC" -> Max[Abs[Lookup[#2, "PctChangeAUC"]]]|> &,
    grouped
];
sensitivitySorted = ReverseSortBy[sensitivity, #MaxAbsPctChangeAUC &];

rankPath = FileNameJoin[{"Results", "experiment7_sensitivity_ranking.csv"}];
rankColumnNames = Keys[sensitivitySorted[[1]]];
rankCsvRows = Prepend[Values /@ sensitivitySorted, rankColumnNames];
Export[rankPath, rankCsvRows, "CSV"];
Print["Saved: ", rankPath];

Print["Top 5 most sensitive parameters: "];
Do[Print["  ", s["Parameter"], ": ", s["MaxAbsPctChangeAUC"], "% max change in PestAUC"], {s, sensitivitySorted[[1 ;; 5]]}];

(* Tornado chart of top 10 *)
top10 = sensitivitySorted[[1 ;; Min[10, Length[sensitivitySorted]]]];
tornadoData = Table[
    Module[{name = t["Parameter"], lowRow, highRow},
        lowRow = SelectFirst[results, #Parameter == name && #Factor == 0.9 &];
        highRow = SelectFirst[results, #Parameter == name && #Factor == 1.1 &];
        {name, lowRow["PctChangeAUC"], highRow["PctChangeAUC"]}
    ],
    {t, top10}
];
tornado = BarChart[
    Table[{tornadoData[[k, 2]], tornadoData[[k, 3]]}, {k, Length[tornadoData]}],
    ChartLabels -> {Table[tornadoData[[k, 1]], {k, Length[tornadoData]}], {"-10%", "+10%"}},
    ChartLegends -> {"param \[Times] 0.9", "param \[Times] 1.1"},
    BarSpacing -> {0, 0.6},
    PlotLabel -> "Sensitivity of PestAUC to \[PlusMinus]10% parameter perturbations (top 10, sorted)",
    FrameLabel -> {"Parameter", "% change in PestAUC"},
    Frame -> True, GridLines -> Automatic, GridLinesStyle -> Directive[Dashed],
    ImageSize -> Large
];
pngPath = FileNameJoin[{"Results", "experiment7_tornado.png"}];
Export[pngPath, tornado];
Print["Saved: ", pngPath];
