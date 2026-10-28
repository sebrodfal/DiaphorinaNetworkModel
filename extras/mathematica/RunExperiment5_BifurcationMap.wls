#!/usr/bin/env wolframscript
(* Experiment 5a: bifurcation map of the SINGLE-NODE impulsive model over
   (delta, T). Long-run diagnostics (200+ release cycles) revealed that the
   demo parameters (delta=35, T=15) do NOT converge to a simple period-T
   orbit: the post-jump state sequence s(nT) keeps wandering with gaps of
   40-80 population units even after 200 cycles, and two different, both
   high-accuracy integrators (ExplicitRungeKutta vs adaptive with
   AccuracyGoal/PrecisionGoal->12) diverge from each other exponentially --
   the signature of deterministic chaos, not a numerical artifact.

   This matters for two reasons: (1) it means "sustained suppression" claims
   for the base model are only established over short horizons (~100 time
   units, ~7 cycles) in the existing notebooks; (2) the network Master
   Stability Function analysis (NetworkStability.wl) requires linearizing
   around a genuine period-T orbit, so we need a (delta,T) region where one
   actually exists before that analysis is meaningful.

   This script classifies each (delta,T) combination as "regular" (settles to
   a period-T, or low-order period-kT, cycle) or "irregular" (chaotic/wide
   dispersion) using the dispersion of the last 30 post-jump states out of a
   120-cycle run, and saves a classification CSV + heatmap PNG. *)

SetDirectory[DirectoryName[$InputFileName]];
Get["NetworkStability.wl"];
If[! DirectoryQ["Results"], CreateDirectory["Results"]];

p1 = {-0.20, 0.040, -0.010, 5.55*10^-4, 0.005};
p2 = {-0.72, -0.010, 0.054, 5.55*10^-4, 0.005};
p3 = {-0.30, 0.010, 5.55*10^-4, 0.005};
p4 = {0.72, -0.072, 5.55*10^-4, 0.005};
params = {p1, p2, p3, p4};
f[{y1_, y2_, y3_, y4_}] := LocalField[params][{y1, y2, y3, y4}];

nCycles = 120;
nTail = 30;

ClassifyRegime[deltaAmt_, Tval_] := Module[{tt, eqs, sol, samples, tail, dispersion},
    tt = nCycles*Tval;
    eqs = {
        s1'[t] == f[{s1[t], s2[t], s3[t], s4[t]}][[1]],
        s2'[t] == f[{s1[t], s2[t], s3[t], s4[t]}][[2]],
        s3'[t] == f[{s1[t], s2[t], s3[t], s4[t]}][[3]],
        s4'[t] == f[{s1[t], s2[t], s3[t], s4[t]}][[4]],
        s1[0] == 50, s2[0] == 10, s3[0] == 15, s4[0] == 30,
        WhenEvent[Mod[t, Tval] == 0 && t > 0, s3[t] -> s3[t] + deltaAmt]
    };
    sol = Quiet[First @ NDSolve[eqs, {s1, s2, s3, s4}, {t, 0, tt},
        Method -> {"TimeIntegration" -> "ExplicitRungeKutta"}]];
    samples = Table[
        {s1[n*Tval + 10^-6], s2[n*Tval + 10^-6], s3[n*Tval + 10^-6], s4[n*Tval + 10^-6]} /. sol,
        {n, nCycles - nTail, nCycles - 1}
    ];
    tail = samples;
    dispersion = Mean[Table[Norm[tail[[k]] - Mean[tail]], {k, Length[tail]}]];
    <|"Delta" -> deltaAmt, "T" -> Tval, "Dispersion" -> dispersion,
      "Regime" -> If[dispersion < 1., "Regular", "Irregular"],
      "s0" -> tail[[-1]]|>
];

deltaGrid = {5, 10, 15, 20, 25, 30, 35, 40};
Tgrid = {5, 8, 11, 15, 20, 25, 30, 40};

results = Flatten @ Table[
    ClassifyRegime[d, tp],
    {d, deltaGrid}, {tp, Tgrid}
];

Print["Completed ", Length[results], " classifications."];

csvPath = FileNameJoin[{"Results", "experiment5_bifurcation_map.csv"}];
tableForCsv = Select[Keys[results[[1]]], # =!= "s0" &];
csvRows = Prepend[Values /@ (KeyDrop[#, "s0"] & /@ results), tableForCsv];
Export[csvPath, csvRows, "CSV"];
Print["Saved: ", csvPath];

(* Heatmap of dispersion, log-scaled for visibility, T on x-axis, delta on y-axis *)
dispByIndex = Table[
    SelectFirst[results, #Delta == d && #T == tp &]["Dispersion"],
    {d, deltaGrid}, {tp, Tgrid}
];
logDisp = Log10[dispByIndex + 10^-3];

hm = ArrayPlot[
    Reverse[logDisp],
    ColorFunction -> "TemperatureMap",
    FrameTicks -> {
        {Table[{Length[deltaGrid] - i + 1, deltaGrid[[i]]}, {i, Length[deltaGrid]}], None},
        {Table[{j, Tgrid[[j]]}, {j, Length[Tgrid]}], None}
    },
    FrameLabel -> {"T (release period)", "\[Delta] (release amount)"},
    PlotLabel -> "log10(dispersion of last 30 post-jump states) \[LongDash] low = periodic, high = irregular/chaotic",
    ImageSize -> Large
];
pngPath = FileNameJoin[{"Results", "experiment5_bifurcation_heatmap.png"}];
Export[pngPath, hm];
Print["Saved: ", pngPath];

regularPoints = Select[results, #Regime == "Regular" &];
Print["Regular (delta,T) points: ", {#Delta, #T} & /@ regularPoints];
