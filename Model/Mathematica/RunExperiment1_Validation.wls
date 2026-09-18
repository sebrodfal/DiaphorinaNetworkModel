#!/usr/bin/env wolframscript
(* Sanity check: with Dd=Dt=0, every node of the network model must evolve
   exactly like the single-node model (SystemEvolve in Functions.m), since the
   diffusive coupling is switched off. This validates that SolveNetworkData
   reduces correctly to the base ODE system before it is used for experiments. *)

SetDirectory[DirectoryName[$InputFileName]];
Get["Functions.m"];
Get["NetworkExperiments.wl"];

p1 = {-0.20, 0.040, -0.010, 5.55*10^-4, 0.005};
p2 = {-0.72, -0.010, 0.054, 5.55*10^-4, 0.005};
p3 = {-0.30, 0.010, 5.55*10^-4, 0.005};
p4 = {0.72, -0.072, 5.55*10^-4, 0.005};
params = {p1, p2, p3, p4};

tt = 100;
n = 3;
ics = {ConstantArray[50, n], ConstantArray[10, n], ConstantArray[15, n], ConstantArray[30, n], tt};

(* Reference: single-node solution, no releases (\[Delta]=0) *)
refEqs = {
    x1s'[t] == x1s[t]*(p1[[1]] + p1[[2]]*x2s[t] + p1[[3]]*x3s[t] - (p1[[4]] + p1[[5]]*(p1[[2]]*x2s[t] + p1[[3]]*x3s[t]))*x1s[t]),
    x2s'[t] == x2s[t]*(p2[[1]] + p2[[2]]*x1s[t] + p2[[3]]*x4s[t] - (p2[[4]] + p2[[5]]*(p2[[2]]*x1s[t] + p2[[3]]*x4s[t]))*x2s[t]),
    x3s'[t] == x3s[t]*(p3[[1]] + p3[[2]]*x1s[t] - (p3[[3]] + p3[[4]]*p3[[2]]*x1s[t])*x3s[t]),
    x4s'[t] == x4s[t]*(p4[[1]] + p4[[2]]*x2s[t] - (p4[[3]] + p4[[4]]*p4[[2]]*x2s[t])*x4s[t]),
    x1s[0] == 50, x2s[0] == 10, x3s[0] == 15, x4s[0] == 30
};
refSol = First @ NDSolve[refEqs, {x1s, x2s, x3s, x4s}, {t, 0, tt},
    Method -> {"TimeIntegration" -> "ExplicitRungeKutta"}];

(* Network run: 3 disconnected-by-coupling nodes (Dd=Dt=0) on a Path graph, no releases *)
graph = MakeTopology["Path", n];
netData = SolveNetworkData[graph, 0, 0, 0, 15, params, ics];

tvals = Subdivide[0., tt, 200];
maxErrX1 = Max[Table[
    Max[Abs[(NodeIF[netData, x1, i] /@ tvals) - ((x1s /. refSol) /@ tvals)]],
    {i, n}
]];
maxErrX3 = Max[Table[
    Max[Abs[(NodeIF[netData, x3, i] /@ tvals) - ((x3s /. refSol) /@ tvals)]],
    {i, n}
]];

(* All nodes should also be pairwise identical to each other, since Dd=Dt=0
   and every node shares the same parameters/ICs *)
maxCrossNodeErr = Max[Table[
    Max[Abs[(NodeIF[netData, x1, 1] /@ tvals) - (NodeIF[netData, x1, i] /@ tvals)]],
    {i, 2, n}
]];

report = <|
    "MaxErrorX1_vs_SingleNode" -> maxErrX1,
    "MaxErrorX3_vs_SingleNode" -> maxErrX3,
    "MaxCrossNodeError" -> maxCrossNodeErr,
    "Pass" -> (maxErrX1 < 10^-6 && maxErrX3 < 10^-6 && maxCrossNodeErr < 10^-9)
|>;

Print["Validation report: ", report];
If[report["Pass"], Print["PASS: network model reduces to single-node model when Dd=Dt=0."],
   Print["FAIL: network model does NOT reduce to single-node model. Check SolveNetworkData."]];
