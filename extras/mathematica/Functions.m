(* ::Package:: *)

(*Package["DiaphorinaPackage"]*)


(* ::Section:: *)
(*Differential Equation Model*)


ChartReport[init_,param1_,param2_,param3_,param4_]:=Transpose@Dataset[List@@@Normal[Dataset/@Association[
{"init"->Association@Thread[{"x1[0]","x2[0]","x3[0]","x4[0]","t"}->init]},
{"x1"->Association@Thread[{"r1","b12","b13","a1","c1"}->param1]},{"x2"->Association@Thread[{"r2","b21","b24","a2","c2"}->param2]},{"x3"->Association@Thread[{"r3","b31","a3","c3"}->param3]},{"x4"->Association@Thread[{"r4","b42","a4","c4"}->param4]}]]]


ClearAll[SystemEvolve];
Options[SystemEvolve]={"LabelPosition"->{0.97,0.95}};
SystemEvolve[\[Delta]_,T_,{{r1_,b12_,b13_,a1_,c1_},{r2_,b21_,b24_,a2_,c2_},{r3_,b31_,a3_,c3_},{r4_,b42_,a4_,c4_}},{Di0_,B0_,x30_,x40_,tt_},opts:OptionsPattern[]]:=Module[
{x,eqs,x1,x2,x3,x4},
eqs=NDSolve[{
(x1'[t]==x1[t]*(r1+b12*x2[t]+b13*x3[t]-(a1+c1*(b12*x2[t]+b13*x3[t]))*x1[t])), (*Diaphorina citri*)
(x2'[t]==x2[t]*(r2+b21*x1[t]+b24*x4[t]-(a2+c2*(b21*x1[t]+b24*x4[t]))*x2[t])), (*Plant shoots*)
(x3'[t]==x3[t]*(r3+b31*x1[t]-(a3+c3*b31*x1[t])*x3[t])), (*Tamarixia radiata*)
(x4'[t]==x4[t]*(r4+b42*x2[t]-(a4+c4*b42*x2[t])*x4[t])), (*Plant vigor*)
x1[0]==Di0,
x2[0]==B0,
x3[0]==x30,
x4[0]==x40,
WhenEvent[Mod[t,T]==0 && t>0, x3[t]->x3[t]+\[Delta]]
},
{x1,x2,x3,x4},{t,0,tt},Method->{"TimeIntegration"->"ExplicitRungeKutta"}];
Plot[Evaluate[{x1[t],x2[t],x3[t],x4[t]}/.eqs],{t,0,tt},
	PlotLegends->Placed[{"Diaphorina citri","Plant shoots","Tamarixia radiata","Plant vigor"},{OptionValue["LabelPosition"],{Right,Top}}],
	PlotRange->All, PlotRangePadding->{{Automatic,Automatic},{Automatic,Scaled[0.3]}},
	FrameLabel->{"time (a.u.)", "population"},GridLines->Automatic,GridLinesStyle->Directive[Dashed],Frame->True,ImageSize->Large]
]


(* ::Section:: *)
(*Network version (test)*)


ClearAll[SystemEvolveNetwork];

Options[SystemEvolveNetwork] = {"LabelPosition" -> {0.97, 0.95}};

SystemEvolveNetwork[
    graph_?GraphQ,
    Dd_, Dt_,
    \[Delta]_, T_,
    {{r1_, b12_, b13_, a1_, c1_},
     {r2_, b21_, b24_, a2_, c2_},
     {r3_, b31_, a3_, c3_},
     {r4_, b42_, a4_, c4_}},
    {Di0_, B0_, x30_, x40_, tt_},
    opts : OptionsPattern[]
] := Module[
    {L, n, eqs, ics, vars, sol, i, j},
    
    (* Laplacian de la red *)
    L = Normal @ KirchhoffMatrix[graph];
    n = VertexCount[graph];
    
    (* Ecuaciones locales con dispersión en Diaphorina y Tamarixia *)
    eqs = Flatten @ {
        
        (* Diaphorina *)
        Table[
            x1[i]'[t] == 
                x1[i][t] * (r1 + b12*x2[i][t] + b13*x3[i][t] - 
                    (a1 + c1*(b12*x2[i][t] + b13*x3[i][t]))*x1[i][t]) -
                Dd * Sum[L[[i, j]] x1[j][t], {j, n}],
            {i, n}
        ],
        
        (* Brotes *)
        Table[
            x2[i]'[t] == 
                x2[i][t] * (r2 + b21*x1[i][t] + b24*x4[i][t] - 
                    (a2 + c2*(b21*x1[i][t] + b24*x4[i][t]))*x2[i][t]),
            {i, n}
        ],
        
        (* Tamarixia *)
        Table[
            x3[i]'[t] == 
                x3[i][t] * (r3 + b31*x1[i][t] - 
                    (a3 + c3*b31*x1[i][t])*x3[i][t]) -
                Dt * Sum[L[[i, j]] x3[j][t], {j, n}],
            {i, n}
        ],
        
        (* Estímulos *)
        Table[
            x4[i]'[t] == 
                x4[i][t] * (r4 + b42*x2[i][t] - 
                    (a4 + c4*b42*x2[i][t])*x4[i][t]),
            {i, n}
        ]
    };
    
    (* Condiciones iniciales *)
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
    
    sol = NDSolve[
        Join[
            eqs,
            ics,
            {
                WhenEvent[
                    Mod[t, T] == 0 && t > 0,
                    Evaluate @ Table[x3[i][t] -> x3[i][t] + \[Delta], {i, n}]
                ]
            }
        ],
        vars,
        {t, 0, tt},
        Method -> {"TimeIntegration" -> "ExplicitRungeKutta"}
    ];
    
    Table[
        Plot[
            Evaluate[{x1[i][t], x2[i][t], x3[i][t], x4[i][t]} /. sol],
            {t, 0, tt},
            PlotLegends -> Placed[{"Diaphorina citri", "Plant shoots", "Tamarixia radiata", "Plant vigor"},
                {OptionValue["LabelPosition"],{Right,Top}}],
            PlotRange -> All,
            PlotRangePadding -> {{Automatic,Automatic},{Automatic,Scaled[0.3]}},
            FrameLabel -> {"time (a.u.)", "population"},
            GridLines -> Automatic,
            GridLinesStyle -> Directive[Dashed],
            Frame -> True,
            ImageSize -> Large,
            PlotLabel -> "Node #" <> ToString[i]
        ],
        {i, n}
    ]
];


(* ::Section:: *)
(*Core Solver (non-plotting, shared by the review-response analyses below)*)


ClearAll[SolveSystem];
SolveSystem[\[Delta]_,T_,{{r1_,b12_,b13_,a1_,c1_},{r2_,b21_,b24_,a2_,c2_},{r3_,b31_,a3_,c3_},{r4_,b42_,a4_,c4_}},{Di0_,B0_,x30_,x40_,tt_}]:=
NDSolve[{
 x1'[t]==x1[t]*(r1+b12*x2[t]+b13*x3[t]-(a1+c1*(b12*x2[t]+b13*x3[t]))*x1[t]),
 x2'[t]==x2[t]*(r2+b21*x1[t]+b24*x4[t]-(a2+c2*(b21*x1[t]+b24*x4[t]))*x2[t]),
 x3'[t]==x3[t]*(r3+b31*x1[t]-(a3+c3*b31*x1[t])*x3[t]),
 x4'[t]==x4[t]*(r4+b42*x2[t]-(a4+c4*b42*x2[t])*x4[t]),
 x1[0]==Di0, x2[0]==B0, x3[0]==x30, x4[0]==x40,
 WhenEvent[Mod[t,T]==0 && t>0, x3[t]->x3[t]+\[Delta]]
},{x1,x2,x3,x4},{t,0,tt},Method->{"TimeIntegration"->"ExplicitRungeKutta"}]


(* ::Section:: *)
(*Release Schedule Sensitivity: (T,\[Delta]) sweep -- responds to R1.3 / R1.14(partial) / R2.7 / R2.8*)


ClearAll[ReleaseSweep];
Options[ReleaseSweep] = {
  "Metrics"->{"MeanX1","PestDays","FractionBelowThreshold"},
  "Window"->Automatic,"Threshold"->10,"Samples"->300
};
ReleaseSweep[Trange_List,\[Delta]range_List,params_,init_,opts:OptionsPattern[]]:=Module[
{tt=Last[init],win,nS,thr,metricsWanted,perPoint},
win = OptionValue["Window"];
If[win===Automatic, win = {tt/2, tt}];
nS = OptionValue["Samples"];
thr = OptionValue["Threshold"];
metricsWanted = OptionValue["Metrics"];
(* Solve once per (T,\[Delta]) grid point and compute all three metrics from that one
   trajectory -- solving the ODE is the expensive step, so this is both cheaper and
   safer than calling ReleaseSweep three times with a different "Metric" each time. *)
perPoint = Flatten[Table[
  Module[{sol,x1f,vals,valsFull,metricsHere},
    sol = Quiet@Check[SolveSystem[\[Delta],T,params,init],$Failed];
    metricsHere = If[sol===$Failed||sol==={},
      <|"MeanX1"->Missing["NotAvailable"],"PestDays"->Missing["NotAvailable"],
        "FractionBelowThreshold"->Missing["NotAvailable"]|>,
      x1f = x1/.First[sol];
      vals = x1f/@Subdivide[win[[1]],win[[2]],nS];
      valsFull = x1f/@Subdivide[0,tt,nS];
      <|
        "MeanX1"->Mean[vals],
        "PestDays"->Mean[valsFull]*tt,
        "FractionBelowThreshold"->N[Count[vals,v_/;v<thr]/Length[vals]]
      |>
    ];
    {T,\[Delta],metricsHere}
  ],{\[Delta],\[Delta]range},{T,Trange}
],1];
Association@Table[
  metric -> Map[{#[[1]],#[[2]],#[[3]][metric]}&, perPoint],
  {metric,metricsWanted}
]
]

ReleaseSweepMetricLabels = <|
  "MeanX1" -> "Mean Diaphorina citri density",
  "PestDays" -> "Pest-days (\[Integral] x1 dt)",
  "FractionBelowThreshold" -> "Fraction of time below threshold"
|>;

ClearAll[ReleaseSweepPlot];
Options[ReleaseSweepPlot] = {
  "LogColor"->{}, (* list of metric keys, e.g. {"MeanX1"}, whose color axis should show log10(value) --
                      use when the metric spans many orders of magnitude and a linear color scale would
                      wash out all but the largest values *)
  "Highlight"->None, (* a single {T,\[Delta]} pair to mark with a star, e.g. the schedule used elsewhere
                         in the paper, so the reader can see where it falls on the surface *)
  "ShowSamples"->True, (* overlay the actual sampled grid points as small dots, so the reader can tell
                           real data from interpolation *)
  "Boundary"->None (* a list of {T,\[Delta],\[Lambda]} triples (e.g. from InvasionExponentSweep) -- overlays
                       its zero-contour (the eradication boundary) on top of the metric heatmap *)
};
(* Plots exactly the metrics present in `result` -- whatever ReleaseSweep computed via its
   "Metrics" option is what gets plotted here, with no separate "which metric" input to keep
   in sync by hand. Returns an Association with the same keys as `result`, one plot each.
   No PlotLabel: as in SystemEvolve, the figure caption belongs in the LaTeX text, not on the plot. *)
ReleaseSweepPlot[result_Association,opts:OptionsPattern[]]:=Module[{makePlot,logMetrics,highlight,showSamples,boundary},
logMetrics = OptionValue["LogColor"];
highlight = OptionValue["Highlight"];
showSamples = OptionValue["ShowSamples"];
boundary = OptionValue["Boundary"];
makePlot[metric_,data_]:=Module[{clean,useLog,plotData,colorFn,legendFn,epi,bndClean},
  clean = Select[data,NumericQ[#[[3]]]&];
  useLog = MemberQ[logMetrics,metric];
  plotData = If[useLog, {#[[1]],#[[2]],Log10[#[[3]]]}& /@ clean, clean];
  colorFn = "TemperatureMap";
  legendFn = If[useLog,
    BarLegend[{"TemperatureMap",MinMax[plotData[[All,3]]]},
      LegendLabel->Placed["log\!\(\*SubscriptBox[\(10\),\(\\\ \)]\)(value)",Top]],
    Automatic];
  epi = {};
  If[TrueQ[showSamples], AppendTo[epi, {PointSize[0.007],GrayLevel[0.25],Point[Most/@clean]}]];
  If[boundary=!=None,
    bndClean = Select[boundary,NumericQ[#[[3]]]&];
    If[Length[bndClean]>=3,
      AppendTo[epi, ListContourPlot[bndClean,Contours->{0},ContourShading->None,
        ContourStyle->{Black,Thickness[0.006]}][[1]]]
    ]
  ];
  If[highlight=!=None,
    AppendTo[epi, {Text[Style["\[FivePointedStar]",22,Black],highlight]}]
  ];
  ListDensityPlot[plotData,
    ColorFunction->colorFn, PlotLegends->legendFn, InterpolationOrder->1,
    PlotRange->All, (* without this, ListDensityPlot's automatic outlier-trimming heuristic can
                        clip the color range and leave the most extreme grid point (here, the
                        smallest log10 value) unfilled/blank instead of colored *)
    Epilog->epi,
    FrameLabel->{"Release period T","Release magnitude \[Delta]"},
    ImageSize->600]
];
Association@KeyValueMap[#1->makePlot[#1,#2]&, result]
]


(* ::Section:: *)
(*Scenario Comparison Metrics (no control / single introduction / periodic release) -- responds to R1.5*)


ClearAll[ScenarioTrajectories];
(* scenarios: list of {name,\[Delta],T} or {name,\[Delta],T,x30override}.
   \[Delta],T define an optional periodic Tamarixia release (WhenEvent every T, +\[Delta]);
   x30override, if present, replaces init's x3(0) -- e.g. 0 for "no parasitoid" vs
   15 for "single introduction"/"periodic release", matching the paper's own
   Figs. 1-3 (norelease.png has zero Tamarixia; noperiodic.png starts from an
   established population and is never topped up again). This is what lets
   ScenarioComparisonPlot reproduce the paper's scenario figure directly, instead
   of a one-off script bypassing these functions. *)
ScenarioTrajectories[params_,init_,scenarios_List]:=Association@Table[
  With[{name=scn[[1]], \[Delta]=scn[[2]], T=scn[[3]],
        thisInit=If[Length[scn]>=4, ReplacePart[init,3->scn[[4]]], init]},
    name -> (x1/.First[SolveSystem[\[Delta],T,params,thisInit]])
  ],
  {scn,scenarios}
]

ClearAll[ScenarioMetrics];
Options[ScenarioMetrics] = {"RecoveryFraction"->0.9,"Samples"->1000,"Window"->Automatic};
(* "MeanDensity"/"PctReductionMeanDensity" are computed over "Window" (default: second half of
   the horizon, {tt/2,tt} -- matching ReleaseSweep's own default and the paper's stated "mean
   over t in [50,100]"), while "PestDays" (the integral proxy meanFull*tt) and "RecoveryTime"
   use the FULL horizon, since recovery can in principle happen before the reporting window
   starts and pest-days is defined as an integral over all of [0,tt]. Earlier versions used one
   Subdivide[0,tt,...] grid for everything, which made "MeanDensity" silently a full-horizon mean
   instead of the windowed one the table caption claims -- this splits the two, as ReleaseSweep
   already did. *)
ScenarioMetrics[params_,init_,scenarios_List,opts:OptionsPattern[]]:=Module[
{tt=Last[init], win, names, sols, refName, refx1f, tsFull, tsWin, results, frac, nS},
nS = OptionValue["Samples"]; frac = OptionValue["RecoveryFraction"];
win = OptionValue["Window"]; If[win===Automatic, win={tt/2,tt}];
names = scenarios[[All,1]];
sols = ScenarioTrajectories[params,init,scenarios];
refName = First[names]; refx1f = sols[refName];
tsFull = Subdivide[0,tt,nS];
tsWin = Subdivide[win[[1]],win[[2]],nS];
results = Table[
  Module[{x1f,valsFull,refValsFull,valsWin,refValsWin,meanDWin,pestDays,recT,\[Delta]here,tRelease,postIdx,postTs,postRatio,dipPos,afterDip,hit},
    x1f = sols[name];
    valsFull = x1f/@tsFull; refValsFull = refx1f/@tsFull;
    valsWin = x1f/@tsWin; refValsWin = refx1f/@tsWin;
    meanDWin = Mean[valsWin]; pestDays = Mean[valsFull]*tt;
    \[Delta]here = scenarios[[Position[names,name][[1,1]],2]];
    (* Where does "after the release" start? For a periodic scenario (\[Delta]>0) that's the
       first pulse, at t=T. For a scenario whose Tamarixia boost is baked into x3(0)
       instead (\[Delta]=0, e.g. "single introduction"), the release already happened at t=0. *)
    tRelease = If[\[Delta]here>0, scenarios[[Position[names,name][[1,1]],3]], 0];
    recT = If[name===refName, 0,
      postIdx = Position[tsFull, t_/;t>tRelease];
      postTs = Extract[tsFull,postIdx];
      postRatio = Extract[valsFull,postIdx]/Extract[refValsFull,postIdx];
      If[Length[postRatio]==0, Missing["NotEnoughData"],
        dipPos = Ordering[postRatio,1][[1]];
        afterDip = Range[dipPos,Length[postRatio]];
        hit = SelectFirst[afterDip, postRatio[[#]]>=frac &];
        If[MissingQ[hit], Missing["NotRecovered"], postTs[[hit]]]
      ]
    ];
    <|"Scenario"->name,"MeanDensity"->meanDWin,"PestDays"->pestDays,"RecoveryTime"->recT,
      "PctReductionMeanDensity"->If[name===refName,0.,N[100*(1-meanDWin/Mean[refValsWin])]]|>
  ],
  {name,names}
];
Dataset[results]
]

ClearAll[ScenarioComparisonPlot];
ScenarioComparisonPlot[params_,init_,scenarios_List]:=Module[{sols,tt=Last[init]},
sols = ScenarioTrajectories[params,init,scenarios];
Plot[Evaluate[Table[sols[name][t],{name,Keys[sols]}]],{t,0,tt},
  PlotLegends->Placed[Keys[sols],{{0.98,0.97},{Right,Top}}],
  PlotRange->All, PlotRangePadding->{{Automatic,Automatic},{Automatic,Scaled[0.3]}},
  FrameLabel->{"time (a.u.)","Diaphorina citri"},
  GridLines->Automatic, GridLinesStyle->Directive[Dashed],
  Frame->True, ImageSize->Large]
]


(* ::Section:: *)
(*Self-Limitation Term Boundedness / Positivity Check -- responds to R2.6*)


ClearAll[SelfLimitationCheck];
Options[SelfLimitationCheck] = {"Samples"->2000};
SelfLimitationCheck[\[Delta]_,T_,params_,init_,opts:OptionsPattern[]]:=Module[
{tt=Last[init],sol,x1f,x2f,x3f,x4f,ts,
 r1,b12,b13,a1,c1, r2,b21,b24,a2,c2, r3,b31,a3,c3, r4,b42,a4,c4,
 coef1,coef2,coef3,coef4},
{{r1,b12,b13,a1,c1},{r2,b21,b24,a2,c2},{r3,b31,a3,c3},{r4,b42,a4,c4}} = params;
sol = SolveSystem[\[Delta],T,params,init];
{x1f,x2f,x3f,x4f} = {x1,x2,x3,x4}/.First[sol];
ts = Subdivide[0,tt,OptionValue["Samples"]];
coef1 = a1 + c1*(b12*x2f[#]+b13*x3f[#]) & /@ ts;
coef2 = a2 + c2*(b21*x1f[#]+b24*x4f[#]) & /@ ts;
coef3 = a3 + c3*b31*x1f[#] & /@ ts;
coef4 = a4 + c4*b42*x2f[#] & /@ ts;
<|
 "x1"-><|"MinCoefficient"->Min[coef1],"ArgMinTime"->ts[[Ordering[coef1,1][[1]]]],"StaysPositive"->Min[coef1]>0,"FractionNegative"->N[Count[coef1,v_/;v<0]/Length[coef1]]|>,
 "x2"-><|"MinCoefficient"->Min[coef2],"ArgMinTime"->ts[[Ordering[coef2,1][[1]]]],"StaysPositive"->Min[coef2]>0,"FractionNegative"->N[Count[coef2,v_/;v<0]/Length[coef2]]|>,
 "x3"-><|"MinCoefficient"->Min[coef3],"ArgMinTime"->ts[[Ordering[coef3,1][[1]]]],"StaysPositive"->Min[coef3]>0,"FractionNegative"->N[Count[coef3,v_/;v<0]/Length[coef3]]|>,
 "x4"-><|"MinCoefficient"->Min[coef4],"ArgMinTime"->ts[[Ordering[coef4,1][[1]]]],"StaysPositive"->Min[coef4]>0,"FractionNegative"->N[Count[coef4,v_/;v<0]/Length[coef4]]|>
|>
]

ClearAll[SelfLimitationPlot];
SelfLimitationPlot[\[Delta]_,T_,params_,init_]:=Module[
{tt=Last[init],sol,x1f,x2f,x3f,x4f,r1,b12,b13,a1,c1,r2,b21,b24,a2,c2,r3,b31,a3,c3,r4,b42,a4,c4},
{{r1,b12,b13,a1,c1},{r2,b21,b24,a2,c2},{r3,b31,a3,c3},{r4,b42,a4,c4}} = params;
sol = SolveSystem[\[Delta],T,params,init];
{x1f,x2f,x3f,x4f} = {x1,x2,x3,x4}/.First[sol];
Plot[Evaluate[{
   a1 + c1*(b12*x2f[t]+b13*x3f[t]),
   a2 + c2*(b21*x1f[t]+b24*x4f[t]),
   a3 + c3*b31*x1f[t],
   a4 + c4*b42*x2f[t]
 }],{t,0,tt},
 PlotLegends->Placed[{"coef. Diaphorina citri","coef. Plant shoots",
   "coef. Tamarixia radiata","coef. Plant vigor"},{{0.98,0.97},{Right,Top}}],
 PlotRange->All, PlotRangePadding->{{Automatic,Automatic},{Automatic,Scaled[0.35]}},
 FrameLabel->{"time (a.u.)",Row[{"self-limitation coefficient  ",Subscript[a,"i"],"+",Subscript[c,"i"]," \[Sigma]\[Beta] ",Subscript[x,"j"]}]},
 GridLines->Automatic, GridLinesStyle->Directive[Dashed],
 Epilog->{Red,Dashed,Line[{{0,0},{tt,0}}]}, Frame->True, ImageSize->Large]
]


(* ::Section:: *)
(*Saturating (Holling type II) Predation Response -- responds to R2.3*)


ClearAll[SolveSystemHolling];
SolveSystemHolling[\[Delta]_,T_,{{r1_,b12_,b13_,a1_,c1_},{r2_,b21_,b24_,a2_,c2_},{r3_,b31_,a3_,c3_},{r4_,b42_,a4_,c4_}},{k13_,k31_},{Di0_,B0_,x30_,x40_,tt_}]:=Module[
{h13,h31},
h13[x3v_]:=b13*x3v/(1+k13*x3v);
h31[x1v_]:=b31*x1v/(1+k31*x1v);
NDSolve[{
 x1'[t]==x1[t]*(r1+b12*x2[t]+h13[x3[t]]-(a1+c1*(b12*x2[t]+h13[x3[t]]))*x1[t]),
 x2'[t]==x2[t]*(r2+b21*x1[t]+b24*x4[t]-(a2+c2*(b21*x1[t]+b24*x4[t]))*x2[t]),
 x3'[t]==x3[t]*(r3+h31[x1[t]]-(a3+c3*h31[x1[t]])*x3[t]),
 x4'[t]==x4[t]*(r4+b42*x2[t]-(a4+c4*b42*x2[t])*x4[t]),
 x1[0]==Di0, x2[0]==B0, x3[0]==x30, x4[0]==x40,
 WhenEvent[Mod[t,T]==0 && t>0, x3[t]->x3[t]+\[Delta]]
},{x1,x2,x3,x4},{t,0,tt},Method->{"TimeIntegration"->"ExplicitRungeKutta"}]
]

ClearAll[HollingComparisonPlot];
HollingComparisonPlot[\[Delta]_,T_,params_,{k13_,k31_},init_]:=Module[
{tt=Last[init],solBase,solHoll,x1b,x1h,x3b,x3h},
solBase = SolveSystem[\[Delta],T,params,init];
solHoll = SolveSystemHolling[\[Delta],T,params,{k13,k31},init];
{x1b,x3b} = {x1,x3}/.First[solBase];
{x1h,x3h} = {x1,x3}/.First[solHoll];
Plot[Evaluate[{x1b[t],x1h[t],x3b[t],x3h[t]}],{t,0,tt},
 PlotLegends->Placed[{
   "Diaphorina citri mass-action (baseline)",
   Row[{"Diaphorina citri Holling II (",Subscript[k,13],"=",k13,", ",Subscript[k,31],"=",k31,")"}],
   "Tamarixia radiata mass-action (baseline)",
   "Tamarixia radiata Holling II"},{{0.98,0.97},{Right,Top}}],
 PlotStyle->{Blue,{Blue,Dashed},Green,{Green,Dashed}},
 PlotRange->All, PlotRangePadding->{{Automatic,Automatic},{Automatic,Scaled[0.35]}},
 FrameLabel->{"time (a.u.)","population"},
 GridLines->Automatic, GridLinesStyle->Directive[Dashed], Frame->True, ImageSize->Large]
]


(* ::Section:: *)
(*Pest-Free Periodic Solution & Invasion Exponent -- responds to R1.14 / R2.6*)
(*Numerical, ergodic-average analogue of a Floquet exponent for the impulsive system*)


ClearAll[SolvePestFreeSubsystem];
SolvePestFreeSubsystem[\[Delta]_,T_,{{r1_,b12_,b13_,a1_,c1_},{r2_,b21_,b24_,a2_,c2_},{r3_,b31_,a3_,c3_},{r4_,b42_,a4_,c4_}},{B0_,x30_,x40_,tt_}]:=
NDSolve[{
 x2'[t]==x2[t]*(r2+b24*x4[t]-(a2+c2*b24*x4[t])*x2[t]),
 x3'[t]==x3[t]*(r3-a3*x3[t]),
 x4'[t]==x4[t]*(r4+b42*x2[t]-(a4+c4*b42*x2[t])*x4[t]),
 x2[0]==B0, x3[0]==x30, x4[0]==x40,
 WhenEvent[Mod[t,T]==0 && t>0, x3[t]->x3[t]+\[Delta]]
},{x2,x3,x4},{t,0,tt},Method->{"TimeIntegration"->"ExplicitRungeKutta"}]

ClearAll[InvasionExponent];
Options[InvasionExponent] = {"TotalTime"->Automatic,"BurnInFraction"->0.3,"Samples"->3000};
InvasionExponent[params_,init3_,\[Delta]_,T_,opts:OptionsPattern[]]:=Module[
{tt,sol,x2f,x3f,r1,b12,b13,ts,integrand,lambda,burn,half1,half2,conv},
tt = OptionValue["TotalTime"];
If[tt===Automatic, tt = Max[400, 60*T]];
burn = OptionValue["BurnInFraction"]*tt;
sol = SolvePestFreeSubsystem[\[Delta],T,params,Append[init3,tt]];
{x2f,x3f} = {x2,x3}/.First[sol];
r1 = params[[1,1]]; b12=params[[1,2]]; b13=params[[1,3]];
ts = Subdivide[burn,tt,OptionValue["Samples"]];
integrand = (r1 + b12*x2f[#] + b13*x3f[#])& /@ ts;
lambda = Mean[integrand];
half1 = Mean[integrand[[1;;Floor[Length[integrand]/2]]]];
half2 = Mean[integrand[[Floor[Length[integrand]/2]+1;;]]];
conv = Abs[half1-half2];
<|"InvasionExponent"->lambda, "PestFreeStable"->(lambda<0), "TotalTime"->tt,
  "BurnIn"->burn, "ConvergenceGap"->conv, "Reliable"->(conv<0.02)|>
]

ClearAll[InvasionExponentSweep];
Options[InvasionExponentSweep] = Options[InvasionExponent];
InvasionExponentSweep[Trange_List,\[Delta]range_List,params_,init3_,opts:OptionsPattern[]]:=Flatten[
Table[
  Module[{res},
    res = Quiet@Check[InvasionExponent[params,init3,\[Delta],T,opts],$Failed];
    If[res===$Failed, {T,\[Delta],Missing["Failed"]}, {T,\[Delta],res["InvasionExponent"]}]
  ],
  {\[Delta],\[Delta]range},{T,Trange}
],1
]

ClearAll[InvasionExponentSweepPlot];
Options[InvasionExponentSweepPlot] = {"Highlight"->None,"ShowSamples"->True};
(* Color scale is always symmetric and centered at zero (Blue=eradication, Red=persistence) --
   maxAbs is set from the data itself so the scale is exactly as wide as the sampled range,
   never wider (which would make the true extremes look artificially mild). *)
InvasionExponentSweepPlot[data_,opts:OptionsPattern[]]:=Module[{clean,maxAbs,dens,bnd,highlight,showSamples,epi},
highlight = OptionValue["Highlight"];
showSamples = OptionValue["ShowSamples"];
clean = Select[data,NumericQ[#[[3]]]&];
maxAbs = Max[Abs[clean[[All,3]]]];
epi = {};
If[TrueQ[showSamples], AppendTo[epi, {PointSize[0.007],GrayLevel[0.25],Point[Most/@clean]}]];
If[highlight=!=None, AppendTo[epi, {Black,Text[Style["\[FivePointedStar]",22,Black],highlight]}]];
dens = ListDensityPlot[clean,
  ColorFunction->(Blend[{Blue,White,Red},(#+maxAbs)/(2 maxAbs)]&), ColorFunctionScaling->False,
  PlotRange->{Full,Full,{-maxAbs,maxAbs}}, Epilog->epi,
  FrameLabel->{"Release period T","Release magnitude \[Delta]"},
  PlotLegends->Automatic, InterpolationOrder->1, ImageSize->600];
bnd = ListContourPlot[clean, Contours->{0}, ContourShading->None, ContourStyle->{Thick,Black}];
Show[dens,bnd]
]


(* ::Section:: *)
(*Robustness to Imperfect Release Synchrony (jitter) -- responds to R1.11*)


ClearAll[SolveSystemJitter];
SolveSystemJitter[\[Delta]_,T_,jitter_,{{r1_,b12_,b13_,a1_,c1_},{r2_,b21_,b24_,a2_,c2_},{r3_,b31_,a3_,c3_},{r4_,b42_,a4_,c4_}},{Di0_,B0_,x30_,x40_,tt_},seed_:1]:=Module[
{nReleases,releaseTimes,events},
SeedRandom[seed];
nReleases = Floor[tt/T];
releaseTimes = Select[
  T*Range[nReleases] + RandomVariate[UniformDistribution[{-jitter,jitter}],nReleases],
  0<#<tt&
];
events = Map[WhenEvent[t==#,x3[t]->x3[t]+\[Delta]]&, releaseTimes];
NDSolve[
  Join[{
    x1'[t]==x1[t]*(r1+b12*x2[t]+b13*x3[t]-(a1+c1*(b12*x2[t]+b13*x3[t]))*x1[t]),
    x2'[t]==x2[t]*(r2+b21*x1[t]+b24*x4[t]-(a2+c2*(b21*x1[t]+b24*x4[t]))*x2[t]),
    x3'[t]==x3[t]*(r3+b31*x1[t]-(a3+c3*b31*x1[t])*x3[t]),
    x4'[t]==x4[t]*(r4+b42*x2[t]-(a4+c4*b42*x2[t])*x4[t]),
    x1[0]==Di0,x2[0]==B0,x3[0]==x30,x4[0]==x40
  }, events],
  {x1,x2,x3,x4},{t,0,tt},Method->{"TimeIntegration"->"ExplicitRungeKutta"}
]
]

ClearAll[JitterRealizations];
(* Shared by JitterRobustnessPlot and JitterRobustnessStats -- runs the nS jittered solves once
   and samples every realization on a common time grid, so the plot and the summary statistics
   are always computed from the exact same runs (no risk of the figure and the reported numbers
   drifting apart, which is what caused the 8-vs-10 mismatch this replaces). *)
JitterRealizations[\[Delta]_,T_,jitter_,params_,init_,nS_,nSamp_]:=Module[{tt=Last[init],ts,valsMat},
ts = Subdivide[0,tt,nSamp];
valsMat = Table[
  Module[{sol,x1f},
    sol = SolveSystemJitter[\[Delta],T,jitter,params,init,seed];
    x1f = x1/.First[sol];
    x1f/@ts
  ],{seed,1,nS}
];
<|"Times"->ts,"Values"->valsMat|>
]

ClearAll[JitterRobustnessStats];
Options[JitterRobustnessStats] = {"Seeds"->50,"Samples"->400,"Window"->Automatic};
(* Mean +/- 1 s.d. of x1 over the given window (default: second half of the horizon, matching the
   convention used by ReleaseSweep/ScenarioMetrics elsewhere), across the jittered realizations. *)
JitterRobustnessStats[\[Delta]_,T_,jitter_,params_,init_,opts:OptionsPattern[]]:=Module[
{tt=Last[init],nS=OptionValue["Seeds"],nSamp=OptionValue["Samples"],win=OptionValue["Window"],
 realiz,ts,valsMat,inWin,windowMeans},
If[win===Automatic, win={tt/2,tt}];
realiz = JitterRealizations[\[Delta],T,jitter,params,init,nS,nSamp];
ts = realiz["Times"]; valsMat = realiz["Values"];
inWin = Flatten[Position[ts, t_/;win[[1]]<=t<=win[[2]]]];
windowMeans = Mean/@valsMat[[All,inWin]];
<|"Seeds"->nS,"Window"->win,"MeanOfMeans"->Mean[windowMeans],"SDOfMeans"->StandardDeviation[windowMeans]|>
]

ClearAll[JitterRobustnessPlot];
Options[JitterRobustnessPlot] = {"Seeds"->50,"Samples"->400};
(* Plots the jittered-release ensemble as a shaded mean +/- 1 s.d. band (not nS individual gray
   lines) against the exact-period trajectory, with an in-figure legend -- easier to read the
   central tendency at a glance, and scales cleanly to nS=50-100 without the plot turning into an
   illegible tangle of lines. *)
JitterRobustnessPlot[\[Delta]_,T_,jitter_,params_,init_,opts:OptionsPattern[]]:=Module[
{tt=Last[init],solExact,x1e,nS=OptionValue["Seeds"],nSamp=OptionValue["Samples"],
 realiz,ts,valsMat,meanVals,sdVals,bandPts,meanPts,exactPts},
solExact = SolveSystem[\[Delta],T,params,init];
x1e = x1/.First[solExact];
realiz = JitterRealizations[\[Delta],T,jitter,params,init,nS,nSamp];
ts = realiz["Times"]; valsMat = realiz["Values"];
meanVals = Mean[valsMat];
sdVals = StandardDeviation[valsMat];
bandPts = Join[Transpose[{ts,meanVals+sdVals}], Reverse[Transpose[{ts,meanVals-sdVals}]]];
meanPts = Transpose[{ts,meanVals}];
exactPts = Table[{t,x1e[t]},{t,ts}];
Legended[
  Graphics[{
    EdgeForm[None], Opacity[0.4], Gray, Polygon[bandPts],
    Gray, Thickness[0.005], Line[meanPts],
    Black, Thickness[0.005], Line[exactPts]
  },
  Frame->True, FrameLabel->{"time (a.u.)",Row[{"Diaphorina citri (",Subscript[x,1],")"}]},
  GridLines->Automatic, GridLinesStyle->Directive[Dashed],
  PlotRange->All, PlotRangePadding->{{Automatic,Automatic},{Automatic,Scaled[0.28]}},
  AspectRatio->1/GoldenRatio, ImageSize->Large],
  Placed[LineLegend[
    {Directive[Gray,Thickness[0.02]],Directive[Black,Thickness[0.02]]},
    {"Jitter: mean \[PlusMinus] 1 s.d. (N="<>ToString[nS]<>", J="<>ToString[jitter]<>")",
     "Exact period (T="<>ToString[T]<>")"}],
    {{0.98,0.97},{Right,Top}}]
]
]
