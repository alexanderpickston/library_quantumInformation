(* ::Package:: *)

BeginPackage["PhotonicQuantumInformation`"]
 

 
 Print[
"Functions for photonic quantum information processing"
]

Begin["`Private`"]
 
(* __Functions__ *)

Fidelity[x_, y_] := Tr[(x^(1/2).y.x^(1/2))^(1/2)];
FidelityPure[x_, y_] := Sqrt[x\[ConjugateTranspose].y.x];

DensityMatrix = #.#\[ConjugateTranspose] &;
Kron[u__] := KroneckerProduct[u];
 
ProjectionMeasurement[x_, y_] := FullSimplify[Abs[x\[ConjugateTranspose].y]^2][[1, 1]];
 
GetAngle[x_] := N@180*FullSimplify@ArcTan[x[[2]]/x[[1]]]/\[Pi];

(* __BasisStates__ *)

h = {{1},{0}};
v = {{0},{1}};
d = 1/Sqrt[2] (h + v);
a = 1/Sqrt[2] (h - v);
r = 1/Sqrt[2] (h + I*v);
l = 1/Sqrt[2] (h - I*v);

hh = Kron[h, h]; hv = Kron[h, v]; hd = Kron[h, d];
ha = Kron[h, a]; hr = Kron[h, r]; hl = Kron[h, l];
vh = Kron[v, h]; vv = Kron[v, v]; vd = Kron[v, d];
va = Kron[v, a]; vr = Kron[v, r]; vl = Kron[v, l];
dh = Kron[d, h]; dv = Kron[d, v]; dd = Kron[d, d];
da = Kron[d, a]; dr = Kron[d, r]; dl = Kron[d, l];
ah = Kron[a, h]; av = Kron[a, v]; ad = Kron[a, v];
aa = Kron[a, a]; ar = Kron[a, r]; al = Kron[a, l];
rh = Kron[r, h]; rv = Kron[r, v]; rd = Kron[r, v];
ra = Kron[r, a]; rr = Kron[r, r]; rl = Kron[r, l];
lh = Kron[l, h]; lv = Kron[l, v]; ld = Kron[l, v];
la = Kron[l, a]; lr = Kron[l, r]; ll = Kron[l, l];

s0 = {{1, 0}, {0, 1}};
sx = {{0, 1}, {1, 0}};
sy = {{0, -I}, {I, 0}};
sz = {{1, 0}, {0, -1}};


(* __Gates__ *)

cNot = {{1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 0, 1}, {0, 0, 1, 0}};
cH = {{1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1/\[Sqrt]2, 1/\[Sqrt]2}, {0, 0, 1/\[Sqrt]2, -(1/\[Sqrt]2)}};
cX = {{1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 0, 1}, {0, 0, 1, 0}};
cY = {{1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 0, -I}, {0, 0, I, 0}};
cZ = {{1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, -1}};
cI = {{1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}};
swap = {{1, 0, 0, 0}, {0, 0, 1, 0}, {0, 1, 0, 0}, {0, 0, 0, 1}};

(* __MoreFucntions__ *)
 
BitFlip[p_, \[Rho]_] := p.cX.\[Rho].cX + (1 - p) \[Rho];
PhaseFlip[p_, \[Rho]_] := p*sz.\[Rho].sz + (1 - p) \[Rho];
BitPhaseFlip[p_, \[Rho]_] := p*sy.\[Rho].sy + (1 - p) \[Rho];
 
DepolControl[p_, \[Rho]_] := (1 - 3/4*p) \[Rho] + p/4*(cX.\[Rho].cX + cY.\[Rho].cY + cZ.\[Rho].cZ);
 
Depol[\[Eta]_, \[Rho]_] := (1 - \[Eta]) \[Rho] + \[Eta] idn/2;
 
Depol1[p_, \[Rho]_] := (1 - 3/4*p) \[Rho] + p/4*(sx.\[Rho].sx + sy.\[Rho].Y + sz.\[Rho].sz);
 
Depol2[p_, \[Rho]_] := p/2*cI + (1 - p) \[Rho];
 
Dephase[{\[Alpha]x_, \[Alpha]y_, \[Alpha]z_}, p_, \[Rho]_] := (1 - p/2) \[Rho] + p/2 (\[Alpha]x sx.\[Rho].sx + \[Alpha]y sy.\[Rho].sy + \[Alpha]z \sz.\[Rho].sz);
 
DepolMemory[\[Eta]_, \[Rho]_] := \[Eta] \[Rho] + (1 - \[Eta]) idn/2;
 
End[] (* correct syntax to end package *)
EndPackage[] (* correct syntax to end package *)
