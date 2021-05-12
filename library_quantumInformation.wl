(* ::Package:: *)

(* ::Text:: *)
(*Basis States*)


h={{1},{0}}; v={{0},{1}}; 
d=1/Sqrt[2] (h+v); a=1/Sqrt[2] (h-v); 
r=1/Sqrt[2] (h+ I*v); l=1/Sqrt[2] (h- I*v);


(* ::Text:: *)
(*Bell States*)


phiplus=1/\[Sqrt]2 (hh+vv); phiplusdm=DensityMatrix[phiplus];
phiminus=1/\[Sqrt]2 (hh-vv); phiminusdm=DensityMatrix[phiminus];
psiplus=1/\[Sqrt]2 (hv+vh); psiplusdm=DensityMatrix[psiplus];
psiminus=1/\[Sqrt]2 (hv-vh); psiminusdm=DensityMatrix[psiminus];


(* ::Text:: *)
(*Density Matrix*)


DensityMatrix=#.#\[ConjugateTranspose]&;

(* normalizedm--Normalize a density matrix *)
normalizedm[dm_]:=dm/Tr[dm]//N//Chop;

(* puredm--Calculate density matrix of a pure state *)
puredm[state_]:=state.state\[ConjugateTranspose]//N//Chop;

fidelity[dm_,idealdm_]:=Block[{idealdmnormalized},
idealdmnormalized=normalizedm[idealdm];
Tr[MatrixPower[MatrixPower[idealdm,1/2].dm.MatrixPower[idealdm,1/2],1/2]]^2//Chop
];

(* purity--Calculate purity of a density matrix *)
purity[dm_]:=Tr[dm.dm]//N//Chop;


(* ::Text:: *)
(*Kronecker Product*)


Kron[u__]:=KroneckerProduct[u];

hh=Kron[h,h]; hv=Kron[h,v]; hd=Kron[h,d];
ha=Kron[h,a]; hr=Kron[h,r]; hl=Kron[h,l]; 

vh=Kron[v,h]; vv=Kron[v,v];vd=Kron[v,d];
va=Kron[v,a]; vr=Kron[v,r];vl=Kron[v,l];

dh=Kron[d,h]; dv=Kron[d,v]; dd=Kron[d,d];
da=Kron[d,a]; dr=Kron[d,r]; dl=Kron[d,l];

ah=Kron[a,h]; av=Kron[a,v]; ad=Kron[a,v];
aa=Kron[a,a]; ar=Kron[a,r]; al=Kron[a,l];

rh=Kron[r,h]; rv=Kron[r,v]; rd=Kron[r,v];
ra=Kron[r,a]; rr=Kron[r,r]; rl=Kron[r,l];

lh=Kron[l,h]; lv=Kron[l,v]; ld=Kron[l,v];
la=Kron[l,a]; lr=Kron[l,r]; ll=Kron[l,l];


(* ::Text:: *)
(*Pauli Matrices*)


s0={{1,0},{0,1}};
sx={{0,1},{1,0}}; 
sy={{0,-I},{I,0}}; 
sz={{1,0},{0,-1}};


(* ::Text:: *)
(*Gate Operations*)


(* controlled not gate *)
cnot={{1,0,0,0},{0,1,0,0},{0,0,0,1},{0,0,1,0}}; 
\!\(\*
TagBox[
RowBox[{"(", GridBox[{
{"1", "0", "0", "0"},
{"0", "1", "0", "0"},
{"0", "0", "0", "1"},
{"0", "0", "1", "0"}
},
GridBoxAlignment->{"Columns" -> {{Center}}, "ColumnsIndexed" -> {}, "Rows" -> {{Baseline}}, "RowsIndexed" -> {}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.7]}, Offset[0.27999999999999997`]}, "ColumnsIndexed" -> {}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}, "RowsIndexed" -> {}}], ")"}],
Function[BoxForm`e$, MatrixForm[BoxForm`e$]]]\);
(* swap gate *)
swap={{1,0,0,0},{0,0,1,0},{0,1,0,0},{0,0,0,1}};
\!\(\*
TagBox[
RowBox[{"(", GridBox[{
{"1", "0", "0", "0"},
{"0", "0", "1", "0"},
{"0", "1", "0", "0"},
{"0", "0", "0", "1"}
},
GridBoxAlignment->{"Columns" -> {{Center}}, "ColumnsIndexed" -> {}, "Rows" -> {{Baseline}}, "RowsIndexed" -> {}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.7]}, Offset[0.27999999999999997`]}, "ColumnsIndexed" -> {}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}, "RowsIndexed" -> {}}], ")"}],
Function[BoxForm`e$, MatrixForm[BoxForm`e$]]]\);
(* root swap *)
rootswap={{1,0,0,0},{0,1/2*(1+I),1/2*(1-I),0},{0,1/2*(1-I),1/2*(1+I),0},{0,0,0,1}};
\!\(\*
TagBox[
RowBox[{"(", GridBox[{
{"1", "0", "0", "0"},
{"0", 
RowBox[{
FractionBox["1", "2"], "+", 
FractionBox["I", "2"]}], 
RowBox[{
FractionBox["1", "2"], "-", 
FractionBox["I", "2"]}], "0"},
{"0", 
RowBox[{
FractionBox["1", "2"], "-", 
FractionBox["I", "2"]}], 
RowBox[{
FractionBox["1", "2"], "+", 
FractionBox["I", "2"]}], "0"},
{"0", "0", "0", "1"}
},
GridBoxAlignment->{"Columns" -> {{Center}}, "ColumnsIndexed" -> {}, "Rows" -> {{Baseline}}, "RowsIndexed" -> {}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.7]}, Offset[0.27999999999999997`]}, "ColumnsIndexed" -> {}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}, "RowsIndexed" -> {}}], ")"}],
Function[BoxForm`e$, MatrixForm[BoxForm`e$]]]\);
(* hadamard *)
H=1/\[Sqrt]2 {{1,1},{1,-1}}; 
\!\(\*
TagBox[
RowBox[{"(", "", GridBox[{
{
FractionBox["1", 
SqrtBox["2"]], 
FractionBox["1", 
SqrtBox["2"]]},
{
FractionBox["1", 
SqrtBox["2"]], 
RowBox[{"-", 
FractionBox["1", 
SqrtBox["2"]]}]}
},
GridBoxAlignment->{"Columns" -> {{Center}}, "Rows" -> {{Baseline}}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.7]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}}], "", ")"}],
Function[BoxForm`e$, MatrixForm[BoxForm`e$]]]\);
(* controlled hadamard *)
cH={{1,0,0,0},{0,1,0,0},{0,0,1/\[Sqrt]2,1/\[Sqrt]2},{0,0,1/\[Sqrt]2,-(1/\[Sqrt]2)}};
\!\(\*
TagBox[
RowBox[{"(", GridBox[{
{"1", "0", "0", "0"},
{"0", "1", "0", "0"},
{"0", "0", 
FractionBox["1", 
SqrtBox["2"]], 
FractionBox["1", 
SqrtBox["2"]]},
{"0", "0", 
FractionBox["1", 
SqrtBox["2"]], 
RowBox[{"-", 
FractionBox["1", 
SqrtBox["2"]]}]}
},
GridBoxAlignment->{"Columns" -> {{Center}}, "ColumnsIndexed" -> {}, "Rows" -> {{Baseline}}, "RowsIndexed" -> {}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.7]}, Offset[0.27999999999999997`]}, "ColumnsIndexed" -> {}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}, "RowsIndexed" -> {}}], ")"}],
Function[BoxForm`e$, MatrixForm[BoxForm`e$]]]\);
(* controlled X gate *)
csx={{1,0,0,0},{0,1,0,0},{0,0,0,1},{0,0,1,0}};
\!\(\*
TagBox[
RowBox[{"(", GridBox[{
{"1", "0", "0", "0"},
{"0", "1", "0", "0"},
{"0", "0", "0", "1"},
{"0", "0", "1", "0"}
},
GridBoxAlignment->{"Columns" -> {{Center}}, "ColumnsIndexed" -> {}, "Rows" -> {{Baseline}}, "RowsIndexed" -> {}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.7]}, Offset[0.27999999999999997`]}, "ColumnsIndexed" -> {}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}, "RowsIndexed" -> {}}], ")"}],
Function[BoxForm`e$, MatrixForm[BoxForm`e$]]]\);
(* controlled Y gate *)
csy={{1,0,0,0},{0,1,0,0},{0,0,0,-I},{0,0,I,0}};
\!\(\*
TagBox[
RowBox[{"(", GridBox[{
{"1", "0", "0", "0"},
{"0", "1", "0", "0"},
{"0", "0", "0", 
RowBox[{"-", "I"}]},
{"0", "0", "I", "0"}
},
GridBoxAlignment->{"Columns" -> {{Center}}, "ColumnsIndexed" -> {}, "Rows" -> {{Baseline}}, "RowsIndexed" -> {}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.7]}, Offset[0.27999999999999997`]}, "ColumnsIndexed" -> {}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}, "RowsIndexed" -> {}}], ")"}],
Function[BoxForm`e$, MatrixForm[BoxForm`e$]]]\);
(* controlled Z gate *)
csz={{1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,-1}};
\!\(\*
TagBox[
RowBox[{"(", GridBox[{
{"1", "0", "0", "0"},
{"0", "1", "0", "0"},
{"0", "0", "1", "0"},
{"0", "0", "0", 
RowBox[{"-", "1"}]}
},
GridBoxAlignment->{"Columns" -> {{Center}}, "ColumnsIndexed" -> {}, "Rows" -> {{Baseline}}, "RowsIndexed" -> {}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.7]}, Offset[0.27999999999999997`]}, "ColumnsIndexed" -> {}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}, "RowsIndexed" -> {}}], ")"}],
Function[BoxForm`e$, MatrixForm[BoxForm`e$]]]\);
(* four x four identity gate *)
cidn={{1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1}};


(* ::Text:: *)
(*Rotation Operations*)


(*https://en.wikipedia.org/wiki/Jones_calculus*)


qwp[t_]:=1/Sqrt[2] ({
 {1+I Cos[2t], I Sin[2t]},
 {I Sin[2t], 1-I Cos[2t]}
});
hwp[t_]:=Exp[(I \[Pi])/2]({
 {Cos[2t], Sin[2t]},
 {Sin[2t], -Cos[2t]}
});

(* rx--Single qubit X-Rotation *)
rx[\[Theta]_]:={{Cos[\[Theta]/2],-I*Sin[\[Theta]/2]},{-I*Sin[\[Theta]/2],Cos[\[Theta]/2]}};

(* ry--Single qubit Y-Rotation *)
ry[\[Theta]_]:={{Cos[\[Theta]/2],-Sin[\[Theta]/2]},{Sin[\[Theta]/2],Cos[\[Theta]/2]}};

(* rz--Single qubit Z-Rotation *)
rz[\[Theta]_]:={{Exp[-I*\[Theta]/2],0},{0,Exp[I*\[Theta]/2]}};


(* ::Text:: *)
(*Channels*)


(* bit flip *) 
bitflip[p_,\[Rho]_]:=p.sx.\[Rho].sx+(1-p)\[Rho];

(* phase flip *)
phaseflip[p_,\[Rho]_]:=p*sz.\[Rho].sz+(1-p)\[Rho];

(* bit phase flip *)
bitphaseflip[p_,\[Rho]_]:=p*sy.\[Rho].sy+(1-p)\[Rho];

(* contorl depolarising channel *)
depolControl[p_,\[Rho]_]:=(1-3/4*p)\[Rho]+p/4*(csx.\[Rho].csx+csy.\[Rho].csy+csz.\[Rho].csz);

(* depolarising channel - as definded by Nielson *)
depol1[p_,\[Rho]_]:=(1-3/4*p)\[Rho]+p/4*(sx.\[Rho].sx+sy.\[Rho].sy+sz.\[Rho].sz);

(* depolarising channel *)
depol2[p_,\[Rho]_]:=p/2*cidn+(1-p)\[Rho];

(* dephasing channel *)
dephase[{\[Alpha]x_,\[Alpha]y_,\[Alpha]z_},p_,\[Rho]_]:=(1-p/2)\[Rho]+p/2 (\[Alpha]x sx.\[Rho].sx+\[Alpha]y sy.\[Rho].sy+\[Alpha]z sz.\[Rho].sz);
depol[\[Eta]_,\[Rho]_]:=(1-\[Eta]) \[Rho]+\[Eta] idn/2;
depolMemory[\[Eta]_,\[Rho]_]:=\[Eta] \[Rho]+(1-\[Eta]) idn/2;


(* ::Text:: *)
(*Hong-Ou-Mandel Interference functions*)


pHOM[\[Tau]_,\[Tau]0_,\[Sigma]_,cc_,V_]:=cc*(1-V* E^(-(1/2) \[Sigma]^2 (\[Tau]-\[Tau]0)^2) );
ppKTPHOM[\[Tau]_,\[Tau]0_,\[Sigma]_,cc_,V_]:=cc*(1-V*UnitTriangle[(\[Tau]-\[Tau]0)/\[Sigma]]);


(* ::Text:: *)
(*Polarisation HOM*)


(* Insert note as to how data should be imported. Need to also fix to be used with 12 by removing ErrorBar plots and introducing the Around function. *)

(*{{0,-152.55`,163471},{1,-132.55`,63567},{2,-112.55000000000001`,58582},{3,-92.55000000000001`,154899},{4,-72.55000000000001`,192629},{5,-52.55000000000001`,110147},{6,-32.55000000000001`,44648}}*)


PolHomFit[data_]:=Module[{fit,fitSingleBands3s,fitMeanBands3s,paras,plots},
fit=NonlinearModelFit[data[[All,{2,3}]],A*(1-V(Cos[4(u0-u)\[Degree]])),{{A,Max@data[[;;,3]]},{V,1},{u0,0}},u];
fitSingleBands3s[u_]=fit["SinglePredictionBands",ConfidenceLevel->0.9973002039367398`];
fitMeanBands3s[u_]=fit["MeanPredictionBands",ConfidenceLevel->0.9973002039367398`];
paras=fit["ParameterTable"];
plots=Show[ErrorListPlot[{#,ErrorBar[Sqrt[#[[2]]]]}&/@data[[All,{2,3}]],Frame->True,LabelStyle->Directive[Black,FontSize->12]],Plot[{fit[u],fitMeanBands3s[u]},{u,data[[1,2]],data[[-1,2]]},PlotStyle->{{Opacity[1],Blue},{Opacity[0.5],Purple},{Opacity[0.5],Blue}},Filling->{2->{1},3->{1}}],ImageSize->Medium];
{fit,paras,plots}];


(* ::Subtitle:: *)
(*Functions*)


(* ::Text:: *)
(*Projection Measurement*)


ProjectionMeasurement[x_,y_]:=FullSimplify[Abs[x\[ConjugateTranspose].y]^2][[1,1]];


(* ::Text:: *)
(*Density Matrix Plots*)


(* Need to finalise this plot which is in QuantumInformationAP.nb inside GHZ states *)
