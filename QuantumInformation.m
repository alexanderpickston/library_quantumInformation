 
(* __Functions__ *)

Fidelity[\[Rho]1_,\[Rho]2_]:=(Tr[MatrixPower[(MatrixPower[\[Rho]1,1/2].\[Rho]2.MatrixPower[\[Rho]1,1/2]),1/2]])^2
FidelityPure[x_, y_] := Sqrt[x\[ConjugateTranspose].y.x];
Purity[\[Rho]_]:=Tr[\[Rho].\[Rho]];

DensityMatrix:=#.#\[ConjugateTranspose]&
Kron:=KroneckerProduct[##]&

ProjectionMeasurement[x_, y_] := 
  FullSimplify[Abs[x\[ConjugateTranspose].y]^2][[1, 1]];

GetAngle[x_] := N@180*FullSimplify@ArcTan[x[[2]]/x[[1]]]/\[Pi];

(*__BasisStates__*)

H = {{1},{0}};
V = {{0},{1}};
P = (H+V)/Sqrt[2];
M = (H-V)/Sqrt[2];
R = (H+I*V)/Sqrt[2];
L = (H-I*V)/Sqrt[2];

linBasis={H,V};
diagBasis={P,M};
circBasis={R,L};

s0 = {{1, 0}, {0, 1}};
sx = {{0, 1}, {1, 0}};
sy = {{0, -I}, {I, 0}};
sz = {{1, 0}, {0, -1}}; 

GHZ[nQubit_]:=1/Sqrt[2] (Kron@@ConstantArray[h,nQubit]+Kron@@ConstantArray[v,nQubit])/;nQubit>=2
LinearCluster[nQubit_]:=(Dot@@(KP@@#&/@Permutations[{Cphase}~Join~ConstantArray[IM,nQubit-2]])).KP@@ConstantArray[P,nQubit]/;nQubit>=3

(* __Gates__ *)

cNot = {{1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 0, 1}, {0, 0, 1, 0}};
cH = {{1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1/\[Sqrt]2, 1/\[Sqrt]2}, {0, 0, 1/\[Sqrt]2, -(1/\[Sqrt]2)}};
cX = {{1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 0, 1}, {0, 0, 1, 0}};
cY = {{1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 0, -I}, {0, 0, I, 0}};
cZ = {{1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, -1}};
cI = {{1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}};
swap = {{1, 0, 0, 0}, {0, 0, 1, 0}, {0, 1, 0, 0}, {0, 0, 0, 1}};

(* __MoreFunctions__ *)
 
BitFlip[p_, \[Rho]_] := p.cX.\[Rho].cX + (1 - p) \[Rho];
PhaseFlip[p_, \[Rho]_] := p*sz.\[Rho].sz + (1 - p) \[Rho];
BitPhaseFlip[p_, \[Rho]_] := p*sy.\[Rho].sy + (1 - p) \[Rho];
 
DepolControl[p_, \[Rho]_] := (1 - 3/4*p) \[Rho] + p/4*(cX.\[Rho].cX + cY.\[Rho].cY + cZ.\[Rho].cZ);
Depol[\[Eta]_, \[Rho]_] := (1 - \[Eta]) \[Rho] + \[Eta] idn/2;
Depol1[p_, \[Rho]_] := (1 - 3/4*p) \[Rho] + p/4*(sx.\[Rho].sx + sy.\[Rho].Y + sz.\[Rho].sz);
Depol2[p_, \[Rho]_] := p/2*cI + (1 - p) \[Rho];
DepolMemory[\[Eta]_, \[Rho]_] := \[Eta] \[Rho] + (1 - \[Eta]) idn/2;
 
Dephase[{\[Alpha]x_, \[Alpha]y_, \[Alpha]z_}, p_, \[Rho]_] := (1 - p/2) \[Rho] + p/2 (\[Alpha]x sx.\[Rho].sx + \[Alpha]y sy.\[Rho].sy + \[Alpha]z \sz.\[Rho].sz); 

(* __Rotations__ *)

QWP[t_]:=1/Sqrt[2] ({
 {1+I Cos[2t], I Sin[2t]},
 {I Sin[2t], 1-I Cos[2t]}
});
HWP[t_]:=Exp[(I \[Pi])/2]({
 {Cos[2t], Sin[2t]},
 {Sin[2t], -Cos[2t]}
});

(* __PlottingFunctions__ *)

Block[{t0, shift, lt, st, alpha0},
 lt = 0.015;
 st = 0.01;
 dmticks = 
  Join[{#, #, {lt, 0}} & /@ Range[-1, 1, 1/4], {#, "", {st, 0}} & /@ 
    Range[-1, 1, 1/16]];
 ]


Charting`iBarChart3D[{}];
prot = Unprotect@Charting`iBarChart3D;
DownValues@Charting`iBarChart3D = 
  DownValues@Charting`iBarChart3D /. 
   HoldPattern[
     lhs : _["ChartBaseStyle"] = 
      cd_[If[cond_, then_, else_], rest__]] :> (If[cond, 
      lhs = cd[rest];
      If[! MemberQ[lhs, _EdgeForm], AppendTo[lhs, then]], 
      lhs = cd[else, rest]]);
Protect /@ prot;


DensityMatrixPlot[densitymatrx_] := Show[BarChart3D[Re[densitymatrx]
   , ChartLayout -> "Grid"
   , ChartElementFunction -> 
    ChartElementDataFunction["GradientScaleCube", 
     "ColorScheme" -> "GrayYellowTones"]
   , BarSpacing -> {0.2, .2}
   , BaseStyle -> {FontFamily -> "YuMincho", FontSize -> 12, Black}
   , Method -> {"Canvas" -> None}
   , ChartBaseStyle -> EdgeForm[{Thickness[.002], Opacity[1], Black}]
   , Ticks -> {None, None, dmticks}
   , Boxed -> True, 
   BoxStyle -> {Directive[Thickness[0.002], Opacity[.6], Gray, Dashed]}
   , FaceGrids -> None
   , Axes -> {False, False, True}
   , AxesStyle -> {Directive[Black, Opacity[1]], 
     Directive[Black, Opacity[1]], Directive[Black]}
   , BoxRatios -> {1, 1, 50}
   , ImageSize -> Large
   , AspectRatio -> 1
   , ViewPoint -> {1.3, -2.4, 2.}
   ]
  ]

DensityMatrixPlotFull[densitymatrx_] := 
 GraphicsGrid[{{BarChart3D[Re[densitymatrx]
     , ChartLayout -> "Grid"
     , ChartElementFunction -> 
      ChartElementDataFunction["GradientScaleCube", 
       "ColorScheme" -> "GrayYellowTones"]
     , BarSpacing -> {0, 0}
     , BaseStyle -> {FontFamily -> font, FontSize -> fontSize, Black}
     , Method -> {"Canvas" -> None}
     , ChartBaseStyle -> EdgeForm[{Thickness[.002], Opacity[1], Black}]
     , Ticks -> {None, None, dmticks}
     , Boxed -> True, 
     BoxStyle -> {Directive[Thickness[0.002], Opacity[.6], Gray, 
        Dashed]}
     , FaceGrids -> None
     , Axes -> {False, False, True}
     , AxesStyle -> {Directive[Black, Opacity[1]], 
       Directive[Black, Opacity[1]], Directive[Black]}
     , BoxRatios -> {.1, .1, 1}
     , ImageSize -> Large
     , AspectRatio -> 0.8
     , ViewPoint -> {1.3, -2.4, 2.}
     ],
    
    BarChart3D[Im[densitymatrx]
     , ChartLayout -> "Grid"
     , ChartElementFunction -> 
      ChartElementDataFunction["GradientScaleCube", 
       "ColorScheme" -> "GrayYellowTones"]
     , BarSpacing -> {0, 0}
     , BaseStyle -> {FontFamily -> font, FontSize -> fontSize, Black}
     , Method -> {"Canvas" -> None}
     , ChartBaseStyle -> EdgeForm[{Thickness[.002], Opacity[1], Black}]
     , Ticks -> {None, None, dmticks}
     , Boxed -> True, 
     BoxStyle -> {Directive[Thickness[0.002], Opacity[.6], Gray, 
        Dashed]}
     , FaceGrids -> None
     , Axes -> {False, False, True}
     , AxesStyle -> {Directive[Black, Opacity[1]], 
       Directive[Black, Opacity[1]], Directive[Black]}
     , BoxRatios -> {.1, .1, 1}
     , ImageSize -> Large
     , AspectRatio -> 0.8
     , ViewPoint -> {1.3, -2.4, 2.}
     ]}}]
 
End[] (* correct syntax to end package *)
EndPackage[] (* correct syntax to end package *)
