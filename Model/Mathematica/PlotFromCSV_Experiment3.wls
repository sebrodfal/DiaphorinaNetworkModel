#!/usr/bin/env wolframscript
(* Re-draws the Experiment 3 summary figure straight from the saved CSV
   (decoupled from the heavy simulation run). Draws error bars manually with
   Graphics primitives instead of relying on the legacy ErrorBarPlots package
   or version-specific Around rendering, so it's guaranteed to produce an
   actual rendered plot instead of an unevaluated expression. *)

SetDirectory[DirectoryName[$InputFileName]];

summaryPath = FileNameJoin[{"Results", "experiment3_replicates_summary.csv"}];
summary = Import[summaryPath, "Dataset", HeaderLines -> 1];

allTopologies = {"Path", "Cycle", "Star", "Complete", "Random", "SmallWorld", "ScaleFree"};

slice = Normal[summary[Select[Abs[#Dd - 0.03] < 10^-9 && Abs[#Dt - 0.03] < 10^-9 && #Strategy == "Uniform" &]]];
sliceByTopo = Association[#Topology -> # & /@ slice];

meanVals = Table[sliceByTopo[topo]["PestAUC_Mean"], {topo, allTopologies}];
stdVals = Table[sliceByTopo[topo]["PestAUC_Std"], {topo, allTopologies}];
n = Length[allTopologies];
capW = 0.12;

points = Table[{i, meanVals[[i]]}, {i, n}];
bars = Table[Line[{{i, meanVals[[i]] - stdVals[[i]]}, {i, meanVals[[i]] + stdVals[[i]]}}], {i, n}];
caps = Table[{
    Line[{{i - capW, meanVals[[i]] - stdVals[[i]]}, {i + capW, meanVals[[i]] - stdVals[[i]]}}],
    Line[{{i - capW, meanVals[[i]] + stdVals[[i]]}, {i + capW, meanVals[[i]] + stdVals[[i]]}}]
}, {i, n}];

errorPlot = Show[
    ListPlot[points, PlotStyle -> {PointSize[0.018], Black}],
    Graphics[{Thick, RGBColor[0.85, 0.5, 0.1], bars, caps}],
    Frame -> True,
    FrameTicks -> {{Automatic, None}, {Table[{i, allTopologies[[i]]}, {i, n}], None}},
    FrameLabel -> {"Topology", "Pest AUC (network-wide, Uniform release, Dd=Dt=0.03)"},
    PlotLabel -> "Mean \[PlusMinus] std of Pest AUC (10 seeds for stochastic topologies)",
    GridLines -> Automatic, GridLinesStyle -> Directive[Dashed],
    PlotRange -> All, ImageSize -> Large
];

pngPath = FileNameJoin[{"Results", "experiment3_summary_errorbars.png"}];
Export[pngPath, errorPlot];
Print["Saved: ", pngPath];
Print["Means: ", meanVals];
Print["Stds: ", stdVals];
