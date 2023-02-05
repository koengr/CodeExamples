(* ::Package:: *)

 BeginPackage[ "usefulFunctions`"]
  MainFunction::usage = 
	"MainFunction[ x] computes a simple function."

(* -------------------------- *)



Module[ {i,a,n,j,ham, tau, v, gees},

X = {{0,1},{1,0}};
Y = {{0,-I},{I,0}};
Z = {{1,0},{0,-1}};
Id2 = IdentityMatrix[2];
XXYY = {{0,0,0,0},{0,0,2,0},{0,2,0,0},{0,0,0,0}};
H = {{1,1},{1,-1}}/Sqrt[2];
HadY = ConjugateTranspose[ {{1, 1}, {I, -I}} / Sqrt[2] ];
CNOT = {{1,0,0,0},{0,1,0,0},{0,0,0,1},{0,0,1,0}};
CNOTup = IdentityMatrix[4][[{1,4,3,2}]] ;

S = MatrixPower[Z,1/2];
Sx = MatrixPower[ X, 1/2];
Sy = MatrixPower[ Y, 1/2];

Vec[a_, n_ ] := Table[ KroneckerDelta[a+1,i]
, {i,1,n} ];
uniform[n_] := Normalize[Table[1, {x,0,2^n-1}]];
ghz[n_] := Normalize[ Vec[0,2^n] + Vec[2^n-1,2^n] ];
binary1 = {"0","1"};
binary3 = { "000","001","010","011","100","101","110", "111" };
binary4 = Flatten[Table[  binary1[[i]] <> binary3[[j]],  {i,1,Length[binary1]}, {j,1,Length[binary3]}]];

gi[n_] := Table[ Sqrt[i (n-i)]/2, {i,1,(n-1)}]; (* All the nearest neighbour couplings for any #sites n *)
Hg1[n_] := 2DiagonalMatrix[ gi[n], -1 ] +2 DiagonalMatrix[ gi[n], +1]; 
Hg[n_] := Module[{},
  gees = gi[n];
  Sum[ gees[[i+1]] KroneckerProduct[ IdentityMatrix[2^i] , XXYY, IdentityMatrix[2^(n-2-i) ] ], {i,0,n-2} ]
]
Hx[n_] := Sum[ KroneckerProduct[IdentityMatrix[2^i], X, IdentityMatrix[2^(n-1-i)]],{i,0,n-1}]
Hz[n_] := Sum[ KroneckerProduct[IdentityMatrix[2^i], Z, IdentityMatrix[2^(n-1-i)]],{i,0,n-1}]
Hy[n_] := Sum[ KroneckerProduct[IdentityMatrix[2^i], Y, IdentityMatrix[2^(n-1-i)]],{i,0,n-1}]
qop[op_,site_,nsites_]:=Module[{order},
order=Log[2,Length[op]];KroneckerProduct[IdentityMatrix[2^(site-1)],op,IdentityMatrix[2^(nsites-site-(order-1))]]
]

polarForm=Expand[#/.z_?NumericQ:>Abs[z] Exp[I Arg[z]]]&;
SetOptions[Manipulator,Appearance->"Open"];
SetOptions[Manipulate, SaveDefinitions->True];

solveTimeEvolution[ ham_, tau_ ] :=  Module[{sol, dim, startstate, v},
dim = Length[ham];
startstate = Normalize[ Table[1, {i,1,dim}]];
sol = NDSolve[ {I  v'[time] == ham . v[time] /. {t-> tau}, v[0] == startstate}, v, {time,0,tau} ] ;
v  /. sol // Flatten//First
]

EvenOddPositions[n_] := Module[ {sequence},
sequence = {False,True};
For[i= 1, i<n, i++, sequence= Join[ sequence, Not /@ sequence] ];
{ Position[ sequence, False ] // Flatten, Position[ sequence, True ] // Flatten}
]; 

(* Specific for the transmon system *)
dees[gees_, deltas_, etas_] := Table[
2 gees[[i]]^2 / (deltas[[i]] - deltas[[i+1]] + etas[[i+1]]) + 2 gees[[i]]^2 / (deltas[[i+1]] - deltas[[i]] + etas[[i]]),
{i, 1, Length[gees]}];
cees[ gees_, deltas_, etas_]:= Table[ 
 gees[[i]] gees[[i+1]] / (deltas[[i]]-deltas[[i+1]]+etas[[i+1]]) + gees[[i]] gees[[i+1]] / (deltas[[i+2]]-deltas[[i+1]]+etas[[i+1]]), 
{i,1,Length[gees]-1}];
gitime[n_,t_]:= gi[n] * Pi/2 / t;

Hdetune[n_, deltas_] := Sum[ deltas[[idx]]qop[ Z, idx, n], {idx, n} ];
Htritshift[n_, dees_] := Sum[ dees[[idx]] (qop[Z,idx,n] . qop[Z,idx+1,n] - qop[Z,idx,n] - qop[Z,idx+1,n]), {idx, n-1} ];
Hint[n_, gees_] := Sum[ gees[[idx]] qop[XXYY, idx, n] , {idx, n-1} ];
Hsecondint[n_,cees_] := Sum[ cees[[idx]] (qop[X, idx, n].qop[X,idx+2,n] + qop[Y,idx,n] . qop[Y,idx+2,n]) , {idx, n-2} ];
Hindirecthop[n_, cees_] := Sum[ cees[[idx]] (qop[X, idx, n].qop[Z,idx+1,n].qop[X,idx+2,n] + qop[Y,idx,n] .qop[Z,idx+1,n] qop[Y,idx+2,n]) , {idx, n-2} ];

solveTimeEvolution[  ham_, tau_ ] :=  Module[{sol},
sol = NDSolve[ {I  v'[time] == ham . v[time] /. {t-> tau}, v[0] == uniform[Log[2,Length[ham]]]}, v, {time,0,tau} ] ;
v  /. sol // Flatten//First
]

]
(* -------------------------- *)

  EndPackage[]





