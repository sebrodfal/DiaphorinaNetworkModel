#!/usr/bin/env wolframscript
(* Summary figure for Experiment 6, read straight from the saved CSV. *)

SetDirectory[DirectoryName[$InputFileName]];
data = Import[FileNameJoin[{"Results", "experiment6_dose_frequency.csv"}], "Dataset", HeaderLines -> 1];

deltaGrid = {10, 20, 30, 40, 50};
Tgrid = {5, 10, 15, 20, 30};
colors = ColorData[97] /@ Range[Length[Tgrid]];

(* Panel 1: CostEfficiency vs delta, one line per T, topology=Complete, Uniform *)
slice1 = Normal[data[Select[#Topology == "Complete" && #Strategy == "Uniform" &]]];
lines1 = Table[
    {colors[[k]], Line[
        Table[{d, SelectFirst[slice1, #Delta == d && #T == Tgrid[[k]] &]["CostEfficiency"]}, {d, deltaGrid}]
    ]},
    {k, Length[Tgrid]}
];
plot1 = Graphics[
    {Thickness[0.004], lines1,
     Table[{colors[[k]], PointSize[0.012],
        Point[Table[{d, SelectFirst[slice1, #Delta == d && #T == Tgrid[[k]] &]["CostEfficiency"]}, {d, deltaGrid}]]},
        {k, Length[Tgrid]}]
    },
    Frame -> True, FrameLabel -> {"\[Delta] (release amount)", "Cost efficiency (AUC avoided / parasitoid released)"},
    PlotLabel -> "Complete topology, Uniform release: cost efficiency vs dose, by period T",
    GridLines -> Automatic, GridLinesStyle -> Directive[Dashed],
    PlotRange -> All, ImageSize -> Large,
    Epilog -> {Inset[
        Grid[Prepend[Table[{Graphics[{colors[[k]], Disk[]}, ImageSize -> 10], "T=" <> ToString[Tgrid[[k]]]}, {k, Length[Tgrid]}],
            {"", "Period"}], Alignment -> Left],
        Scaled[{0.85, 0.75}]]}
];

(* Panel 2: PestAUC by topology, best (Uniform) vs best (TargetedAtSource) *)
best = Normal[Import[FileNameJoin[{"Results", "experiment6_best_dose_frequency.csv"}], "Dataset", HeaderLines -> 1]];
topologies = {"Path", "Star", "Complete", "Random"};
uniformAUC = Table[SelectFirst[best, #Topology == t && #Strategy == "Uniform" &]["PestAUC"], {t, topologies}];
targetedAUC = Table[SelectFirst[best, #Topology == t && #Strategy == "TargetedAtSource" &]["PestAUC"], {t, topologies}];

plot2 = BarChart[
    {uniformAUC, targetedAUC},
    ChartLabels -> {topologies, {"Uniform", "TargetedAtSource"}},
    ChartLegends -> {"Uniform", "TargetedAtSource"},
    ChartStyle -> {ColorData[97][1], ColorData[97][2]},
    PlotLabel -> "Best-cost-efficiency PestAUC by topology and strategy",
    FrameLabel -> {"Topology", "Pest AUC"},
    Frame -> True, GridLines -> Automatic, GridLinesStyle -> Directive[Dashed],
    ImageSize -> Large
];

combined = Column[{plot1, plot2}, Spacings -> 2];
pngPath = FileNameJoin[{"Results", "experiment6_summary.png"}];
Export[pngPath, combined];
Print["Saved: ", pngPath];
