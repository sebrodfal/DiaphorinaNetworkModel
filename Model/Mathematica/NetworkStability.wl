(* ::Package:: *)

(* Linear (Floquet) stability analysis of the impulsive periodic orbit for the
   DiaphorinaModel network extension.

   Key fact exploited: with a Uniform release strategy and identical local
   parameters at every node, the synchronized state x_i(t) = s(t) for all i
   (s(t) = the single-node periodic orbit) is an exact solution of the network
   system, because graph-Laplacian rows sum to zero so the diffusion term
   vanishes identically on the synchronized manifold. Linearizing the network
   ODE around that orbit and projecting onto the graph Laplacian's eigenbasis
   decouples the 4n-dimensional variational equation into n independent 4x4
   "master stability" problems, one per Laplacian eigenvalue mu_k:

       Y_k'(t) = (Df(s(t)) - mu_k * D) Y_k(t),   Y_k(0) = Identity,   D = diag(Dd,0,Dt,0)

   mu_1 = 0 reproduces the single-node Floquet problem (the synchronized mode
   itself). mu_2..mu_n are the "transverse" modes: if the spectral radius of
   Y_k(T) is < 1 for all of them, small spatial (desynchronizing) perturbations
   decay and the uniform periodic regime is stable on that topology; if any
   >1, the synchronized state loses stability and the network is expected to
   develop persistent spatial heterogeneity.

   Because the jump map (x3 -> x3 + delta, applied uniformly) is an additive
   constant, its Jacobian is the identity, so it does not enter the monodromy
   -- only the continuous-time flow over one period T does.

   MasterMonodromy/MasterStabilityValue/MasterStabilityCurve below require s(t) to
   be a genuine, converged period-T orbit (via FindPeriodicState) -- exact Floquet
   theory needs that. RunExperiment5_BifurcationMap.wls found that, at the base
   single-node parameters, NO (delta,T) combination tried settles to such an orbit
   within 120 release cycles (every point classified "Irregular", including the
   paper's own baseline delta=35,T=15) -- long runs look aperiodic/chaotic instead.

   MasterLyapunovExponent/MasterLyapunovSpectrum generalize the same synchronized-
   manifold argument to that case. Master Stability Function theory (Pecora-Carroll)
   was originally developed for chaotic synchronized attractors, not just periodic
   ones -- FindPeriodicState's exact monodromy over one period T is the special case,
   not a requirement. Instead of integrating Y_k(T) once from a converged s0, these
   integrate a single perturbation vector v(t) along a long reference trajectory
   s(t) (periodic or not), renormalizing every "SegmentLength" to prevent overflow
   (Benettin's method), and report the resulting finite-time Lyapunov exponent of
   that transverse mode -- lambda_k < 0 still means the mode decays (stable),
   lambda_k > 0 still means it grows (desynchronizes), and lambda_k reduces to the
   ordinary Floquet exponent exactly when s(t) genuinely is period-T. Same
   ConvergenceGap/Reliable diagnostic convention as InvasionExponent in the
   single-node ReproduceResults.m (R1.14): split the measurement window in half and
   compare the two halves' exponent estimates. *)


(* ---------- Local vector field and its symbolic Jacobian ---------- *)

ClearAll[LocalField];
LocalField[{{r1_, b12_, b13_, a1_, c1_}, {r2_, b21_, b24_, a2_, c2_},
    {r3_, b31_, a3_, c3_}, {r4_, b42_, a4_, c4_}}][{y1_, y2_, y3_, y4_}] := {
    y1*(r1 + b12*y2 + b13*y3 - (a1 + c1*(b12*y2 + b13*y3))*y1),
    y2*(r2 + b21*y1 + b24*y4 - (a2 + c2*(b21*y1 + b24*y4))*y2),
    y3*(r3 + b31*y1 - (a3 + c3*b31*y1)*y3),
    y4*(r4 + b42*y2 - (a4 + c4*b42*y2)*y4)
};

(* Fixed template symbols (NOT Module-local) used only to build the symbolic
   Jacobian. Module[{y1,y2,y3,y4}, Function[{y1,y2,y3,y4}, Evaluate[J]]] is a
   classic trap: Module renames its locals, but Function does its own
   independent hygiene-renaming of its parameter list, and the two renamings
   don't line up -- the returned function silently fails to substitute its
   arguments. Keeping the template variables outside any Module sidesteps
   that entirely. *)
jacTemplateVars = {ux1, ux2, ux3, ux4};

ClearAll[LocalJacobianFunction];
LocalJacobianFunction[params_] := Module[{f, J},
    f = LocalField[params][jacTemplateVars];
    J = D[f, {jacTemplateVars}];
    Function[{y1, y2, y3, y4}, Evaluate[J /. Thread[jacTemplateVars -> {y1, y2, y3, y4}]]]
];


(* ---------- Reference periodic orbit (single node, converged) ---------- *)

ClearAll[FindPeriodicState];
FindPeriodicState[params_, \[Delta]_, T_, nPeriods_: 25] := Module[
    {f, eqs, sol, s0, sPrev, tBurn, gap},
    f[{y1_, y2_, y3_, y4_}] := LocalField[params][{y1, y2, y3, y4}];
    tBurn = nPeriods*T;
    eqs = {
        s1'[t] == f[{s1[t], s2[t], s3[t], s4[t]}][[1]],
        s2'[t] == f[{s1[t], s2[t], s3[t], s4[t]}][[2]],
        s3'[t] == f[{s1[t], s2[t], s3[t], s4[t]}][[3]],
        s4'[t] == f[{s1[t], s2[t], s3[t], s4[t]}][[4]],
        s1[0] == 50, s2[0] == 10, s3[0] == 15, s4[0] == 30,
        WhenEvent[Mod[t, T] == 0 && t > 0, s3[t] -> s3[t] + \[Delta]]
    };
    sol = First @ NDSolve[eqs, {s1, s2, s3, s4}, {t, 0, tBurn},
        Method -> {"TimeIntegration" -> "ExplicitRungeKutta"}];
    s0 = {s1[tBurn - T + 10^-6], s2[tBurn - T + 10^-6], s3[tBurn - T + 10^-6], s4[tBurn - T + 10^-6]} /. sol;
    sPrev = {s1[tBurn - 2 T + 10^-6], s2[tBurn - 2 T + 10^-6], s3[tBurn - 2 T + 10^-6], s4[tBurn - 2 T + 10^-6]} /. sol;
    gap = Norm[s0 - sPrev];
    <|"s0" -> s0, "ConvergenceGap" -> gap, "Solution" -> sol, "tBurn" -> tBurn|>
];


(* ---------- Master stability: Y_k(T) for a single scalar mu ---------- *)

ClearAll[MasterMonodromy];
MasterMonodromy[params_, Dd_, Dt_, \[Mu]_, s0_, T_] := Module[
    {f, Jf, Dmat, sVars, yVars, sEqs, sICs, yEqs, yICs, sol, Ymat},
    f[{y1_, y2_, y3_, y4_}] := LocalField[params][{y1, y2, y3, y4}];
    Jf = LocalJacobianFunction[params];
    Dmat = DiagonalMatrix[{Dd, 0, Dt, 0}];

    sEqs = {
        s1'[t] == f[{s1[t], s2[t], s3[t], s4[t]}][[1]],
        s2'[t] == f[{s1[t], s2[t], s3[t], s4[t]}][[2]],
        s3'[t] == f[{s1[t], s2[t], s3[t], s4[t]}][[3]],
        s4'[t] == f[{s1[t], s2[t], s3[t], s4[t]}][[4]]
    };
    sICs = Thread[{s1[0], s2[0], s3[0], s4[0]} == s0];

    yEqs = Flatten @ Table[
        y[i, j]'[t] == Sum[
            (Jf[s1[t], s2[t], s3[t], s4[t]][[i, k]] - \[Mu]*Dmat[[i, k]]) * y[k, j][t],
            {k, 1, 4}
        ],
        {i, 1, 4}, {j, 1, 4}
    ];
    yICs = Flatten @ Table[y[i, j][0] == If[i == j, 1, 0], {i, 1, 4}, {j, 1, 4}];

    sol = First @ NDSolve[
        Join[sEqs, sICs, yEqs, yICs],
        Join[{s1, s2, s3, s4}, Flatten @ Table[y[i, j], {i, 1, 4}, {j, 1, 4}]],
        {t, 0, T},
        Method -> {"TimeIntegration" -> "ExplicitRungeKutta"}
    ];

    Ymat = Table[(y[i, j] /. sol)[T], {i, 1, 4}, {j, 1, 4}]
];

ClearAll[MasterStabilityValue];
MasterStabilityValue[params_, Dd_, Dt_, \[Mu]_, s0_, T_] :=
    Max[Abs[Eigenvalues[MasterMonodromy[params, Dd, Dt, \[Mu], s0, T]]]];

ClearAll[MasterStabilityCurve];
MasterStabilityCurve[params_, Dd_, Dt_, s0_, T_, muValues_List] :=
    Table[<|"mu" -> \[Mu], "SpectralRadius" -> MasterStabilityValue[params, Dd, Dt, \[Mu], s0, T]|>,
        {\[Mu], muValues}];


(* ---------- Master stability without assuming periodicity: finite-time Lyapunov
   exponent of the transverse mode, via Benettin renormalization ---------- *)

ClearAll[ReferenceTrajectory];
(* Solves the local (single-node) system over [0,tTotal] -- this stands in for
   FindPeriodicState's converged s0 without requiring convergence to an exact
   periodic orbit; it is simply the trajectory the whole (synchronized) network
   would follow if every node started at `init` and stayed perfectly in sync. *)
ReferenceTrajectory[params_, \[Delta]_, T_, tTotal_, init_: {50, 10, 15, 30}] := Module[
    {f, events, sol},
    f[{y1_, y2_, y3_, y4_}] := LocalField[params][{y1, y2, y3, y4}];
    events = If[\[Delta] == 0, {}, {WhenEvent[Mod[t, T] == 0 && t > 0, rs3[t] -> rs3[t] + \[Delta]]}];
    sol = First @ NDSolve[
        Join[{
            rs1'[t] == f[{rs1[t], rs2[t], rs3[t], rs4[t]}][[1]],
            rs2'[t] == f[{rs1[t], rs2[t], rs3[t], rs4[t]}][[2]],
            rs3'[t] == f[{rs1[t], rs2[t], rs3[t], rs4[t]}][[3]],
            rs4'[t] == f[{rs1[t], rs2[t], rs3[t], rs4[t]}][[4]],
            rs1[0] == init[[1]], rs2[0] == init[[2]], rs3[0] == init[[3]], rs4[0] == init[[4]]
        }, events],
        {rs1, rs2, rs3, rs4}, {t, 0, tTotal},
        Method -> {"TimeIntegration" -> "ExplicitRungeKutta"}
    ];
    <|"s1" -> rs1 /. sol, "s2" -> rs2 /. sol, "s3" -> rs3 /. sol, "s4" -> rs4 /. sol|>
];

ClearAll[PerturbationSegment];
(* Integrates the transverse variational equation for a single vector v (not the
   full 4x4 fundamental matrix) across one segment [ta,tb], using the already-solved
   reference functions sFns={s1,s2,s3,s4} as time-dependent coefficients. Tracking
   one vector instead of the matrix is the standard simplification for the LARGEST
   Lyapunov exponent (Benettin et al. 1980): a generic v converges to the dominant
   growth direction after a few segments regardless of where it started. *)
PerturbationSegment[Jf_, Dmat_, \[Mu]_, sFns : {_, _, _, _}, {ta_, tb_}, v0 : {_, _, _, _}] := Module[
    {s1f, s2f, s3f, s4f, sol},
    {s1f, s2f, s3f, s4f} = sFns;
    sol = First @ NDSolve[
        Join[
            Table[
                v[k]'[t] == Sum[
                    (Jf[s1f[t], s2f[t], s3f[t], s4f[t]][[k, j]] - \[Mu]*Dmat[[k, j]]) * v[j][t],
                    {j, 1, 4}
                ],
                {k, 1, 4}
            ],
            Table[v[k][ta] == v0[[k]], {k, 1, 4}]
        ],
        Table[v[k], {k, 1, 4}], {t, ta, tb},
        Method -> {"TimeIntegration" -> "ExplicitRungeKutta"}
    ];
    Table[(v[k] /. sol)[tb], {k, 1, 4}]
];

ClearAll[LyapunovExponentFromReference];
Options[LyapunovExponentFromReference] = {"ConvergenceThreshold" -> 0.02};
(* Shared core: given an already-solved reference trajectory (sFns), Benettin-
   renormalizes a single perturbation vector over nSeg segments of length segLen,
   accumulating log-growth functionally via NestList (no explicit loop) instead of
   mutating an accumulator in a Do/While -- each state is {time, unit vector,
   cumulative log-sum}, and the step function is a pure function of the previous
   state only. *)
LyapunovExponentFromReference[Jf_, Dmat_, \[Mu]_, sFns : {_, _, _, _}, burn_, segLen_, nSeg_Integer, v0 : {_, _, _, _}, opts : OptionsPattern[]] := Module[
    {step, states, logSums, \[Lambda], half, \[Lambda]First, \[Lambda]Second, gap},
    step[{tC_, vC_, logSum_}] := Module[{vEnd, nrm},
        vEnd = PerturbationSegment[Jf, Dmat, \[Mu], sFns, {tC, tC + segLen}, vC];
        nrm = Norm[vEnd];
        {tC + segLen, vEnd/nrm, logSum + Log[nrm]}
    ];
    states = NestList[step, {burn, Normalize[v0], 0.}, nSeg];
    logSums = states[[All, 3]];
    \[Lambda] = Last[logSums]/(nSeg*segLen);
    half = Floor[nSeg/2];
    \[Lambda]First = If[half >= 1, logSums[[half + 1]]/(half*segLen), Missing["NotEnoughSegments"]];
    \[Lambda]Second = If[nSeg - half >= 1,
        (Last[logSums] - logSums[[half + 1]])/((nSeg - half)*segLen),
        Missing["NotEnoughSegments"]
    ];
    gap = If[NumericQ[\[Lambda]First] && NumericQ[\[Lambda]Second], Abs[\[Lambda]First - \[Lambda]Second], Missing["NotAvailable"]];
    <|"LyapunovExponent" -> \[Lambda], "ConvergenceGap" -> gap,
      "Reliable" -> (NumericQ[gap] && gap < OptionValue["ConvergenceThreshold"])|>
];

ClearAll[MasterLyapunovExponent];
Options[MasterLyapunovExponent] = {
    "TotalTime" -> Automatic,       (* default: long enough for ~30-40 segments even at large T *)
    "BurnInFraction" -> 0.3,        (* same convention as InvasionExponent *)
    "SegmentLength" -> Automatic,   (* default: one release period T (or 10 if delta=0) *)
    "InitialPerturbation" -> {1., 0., 0., 0.},
    "ConvergenceThreshold" -> 0.02,
    "Init" -> {50, 10, 15, 30}
};
MasterLyapunovExponent[params_, Dd_, Dt_, \[Mu]_, \[Delta]_, T_, opts : OptionsPattern[]] := Module[
    {tTotal, burn, segLen, nSeg, Dmat, Jf, ref, sFns, res},
    tTotal = OptionValue["TotalTime"];
    If[tTotal === Automatic, tTotal = Max[600, 80*If[T > 0, T, 10]]];
    burn = OptionValue["BurnInFraction"]*tTotal;
    segLen = OptionValue["SegmentLength"];
    If[segLen === Automatic, segLen = If[T > 0, T, 10]];
    nSeg = Floor[(tTotal - burn)/segLen];
    Dmat = DiagonalMatrix[{Dd, 0, Dt, 0}];
    Jf = LocalJacobianFunction[params];
    ref = ReferenceTrajectory[params, \[Delta], T, tTotal, OptionValue["Init"]];
    sFns = {ref["s1"], ref["s2"], ref["s3"], ref["s4"]};
    res = LyapunovExponentFromReference[Jf, Dmat, \[Mu], sFns, burn, segLen, nSeg,
        OptionValue["InitialPerturbation"], "ConvergenceThreshold" -> OptionValue["ConvergenceThreshold"]];
    Join[res, <|"ModeStable" -> (res["LyapunovExponent"] < 0), "TotalTime" -> tTotal,
        "BurnIn" -> burn, "NumSegments" -> nSeg, "SegmentLength" -> segLen|>]
];

ClearAll[MasterLyapunovSpectrum];
Options[MasterLyapunovSpectrum] = Options[MasterLyapunovExponent];
(* Solves the reference trajectory ONCE and reuses it for every mu -- the expensive
   step is ReferenceTrajectory, not the (cheap, 4-variable) perturbation segments,
   so this is far cheaper than calling MasterLyapunovExponent once per mu. *)
MasterLyapunovSpectrum[params_, Dd_, Dt_, muValues_List, \[Delta]_, T_, opts : OptionsPattern[]] := Module[
    {tTotal, burn, segLen, nSeg, Dmat, Jf, ref, sFns},
    tTotal = OptionValue["TotalTime"];
    If[tTotal === Automatic, tTotal = Max[600, 80*If[T > 0, T, 10]]];
    burn = OptionValue["BurnInFraction"]*tTotal;
    segLen = OptionValue["SegmentLength"];
    If[segLen === Automatic, segLen = If[T > 0, T, 10]];
    nSeg = Floor[(tTotal - burn)/segLen];
    Dmat = DiagonalMatrix[{Dd, 0, Dt, 0}];
    Jf = LocalJacobianFunction[params];
    ref = ReferenceTrajectory[params, \[Delta], T, tTotal, OptionValue["Init"]];
    sFns = {ref["s1"], ref["s2"], ref["s3"], ref["s4"]};
    Table[
        Module[{res},
            res = LyapunovExponentFromReference[Jf, Dmat, \[Mu], sFns, burn, segLen, nSeg,
                OptionValue["InitialPerturbation"], "ConvergenceThreshold" -> OptionValue["ConvergenceThreshold"]];
            Join[res, <|"mu" -> \[Mu], "ModeStable" -> res["LyapunovExponent"] < 0|>]
        ],
        {\[Mu], muValues}
    ]
];


(* ---------- Brute-force full-network monodromy (finite differences) ---------- *)
(* Ground truth used only to validate MasterMonodromy; O(4n) integrations, so
   keep n small when calling this. *)

ClearAll[ContinuousFlowNetwork];
ContinuousFlowNetwork[graph_?GraphQ, Dd_, Dt_, params_, X0flat_, T_] := Module[
    {L, n, eqs, ics, vars, sol, x1v, x2v, x3v, x4v},
    L = Normal @ KirchhoffMatrix[graph];
    n = VertexCount[graph];
    {x1v, x2v, x3v, x4v} = TakeList[X0flat, {n, n, n, n}];
    eqs = Flatten @ {
        Table[q1[i]'[t] == LocalField[params][{q1[i][t], q2[i][t], q3[i][t], q4[i][t]}][[1]] -
            Dd*Sum[L[[i, j]] q1[j][t], {j, n}], {i, n}],
        Table[q2[i]'[t] == LocalField[params][{q1[i][t], q2[i][t], q3[i][t], q4[i][t]}][[2]], {i, n}],
        Table[q3[i]'[t] == LocalField[params][{q1[i][t], q2[i][t], q3[i][t], q4[i][t]}][[3]] -
            Dt*Sum[L[[i, j]] q3[j][t], {j, n}], {i, n}],
        Table[q4[i]'[t] == LocalField[params][{q1[i][t], q2[i][t], q3[i][t], q4[i][t]}][[4]], {i, n}]
    };
    ics = Flatten @ {
        Table[q1[i][0] == x1v[[i]], {i, n}], Table[q2[i][0] == x2v[[i]], {i, n}],
        Table[q3[i][0] == x3v[[i]], {i, n}], Table[q4[i][0] == x4v[[i]], {i, n}]
    };
    vars = Flatten @ {Table[q1[i], {i, n}], Table[q2[i], {i, n}], Table[q3[i], {i, n}], Table[q4[i], {i, n}]};
    sol = First @ NDSolve[Join[eqs, ics], vars, {t, 0, T},
        Method -> {"TimeIntegration" -> "ExplicitRungeKutta"}];
    Flatten[Table[(v[i] /. sol)[T], {v, {q1, q2, q3, q4}}, {i, n}]]
];

ClearAll[NetworkMonodromy];
NetworkMonodromy[graph_?GraphQ, Dd_, Dt_, params_, X0flat_, T_, eps_: 10^-5] := Module[
    {dim, F0, cols},
    dim = Length[X0flat];
    F0 = ContinuousFlowNetwork[graph, Dd, Dt, params, X0flat, T];
    cols = Table[
        (ContinuousFlowNetwork[graph, Dd, Dt, params, X0flat + eps*UnitVector[dim, d], T] - F0)/eps,
        {d, dim}
    ];
    Transpose[cols]
];

ClearAll[NetworkFloquetMultipliers];
NetworkFloquetMultipliers[graph_?GraphQ, Dd_, Dt_, params_, X0flat_, T_, eps_: 10^-5] :=
    Eigenvalues[NetworkMonodromy[graph, Dd, Dt, params, X0flat, T, eps]];
