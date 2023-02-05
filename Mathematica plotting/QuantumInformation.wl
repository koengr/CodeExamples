(* ::Package:: *)

(* ::Section:: *)
(*Public section: Common definitions*)


(* ::Text:: *)
(*Public section --- Everything defined here is directly evaluated. Any parameters referenced will turn black *)


BeginPackage[ "QuantumInformation`"];

(* Plotting preferences *)
SetOptions[Plot, ImageSize -> {400, 300}, PlotLegends -> Automatic];
SetOptions[ ListPlot, ImageSize -> {400, 300}, Joined -> True, 
  PlotMarkers -> {Automatic, Tiny}, PlotRange -> All, 
  PlotLegends -> Automatic ];
SetOptions[ListLinePlot, ImageSize -> {400, 300}];
SetOptions[Histogram, ImageSize -> {400, 300}];
SetOptions[ListContourPlot, ImageSize -> {400, 300}];
SetOptions[Manipulator, Appearance -> "Open"];

(* Standard quantum information *)
splus = {{0, 0}, {1, 0}};
smin = ConjugateTranspose[splus];
X = {{0, 1}, {1, 0}};
Y = {{0, -I}, {I, 0}};
Z = {{1, 0}, {0, -1}};
H = (X+Z)/Sqrt[2];
Id2 = IdentityMatrix[2];
XXYY = {{0, 0, 0, 0}, {0, 0, 2, 0}, {0, 2, 0, 0}, {0, 0, 0, 0}};

(* Define all functions / modules *)
allstates::usage = "allstates[nlevels, nqubits] gives all nlevels-airy numbers (e.g. binary) on nqubits digits";
particlestates::usage = "particlestates[nlevels, nqubits, nparticles] gives the fock representations corresponding to nparticles";
particlestatestr::usage = "\!\(\*
StyleBox[\"particlestates\",\nFontWeight->\"Bold\"]\)[\!\(\*
StyleBox[\"nlevels\",\nFontSlant->\"Italic\"]\), \!\(\*
StyleBox[\"nqubits\",\nFontSlant->\"Italic\"]\), \!\(\*
StyleBox[\"nparticles\",\nFontSlant->\"Italic\"]\)] gives all states with nparticles as compact strings";
allstatestr::usage = "allstatestr[nlevels, ntransmons] gives a string representing all fock states";
particleindices::usage = "particleindices[nlevels, ntransmons, nparticles] gives the indices of nparticle states";


minslices::usage = "minslices[tmax, w, minsliceperosc:8 ] calculates the minimum number of timeslices to make oscillations with angular frequency w noticable: minsliceperosc * tmax * w / (2 Pi) "

qop; Vec; matrixBasis; 
evolveState; evolveState2;
evolveExp; evolveExpFloquet;
evolveUnitary; evolveUnitary2;
 interpolationElement; makePlotGrid; pidx;



(* ::Section:: *)
(*Private section: All more elaborate functions go here*)


(* Private section. Any defined functions/modules should placed here *)
Begin["`Private`"];

(*  Get a list of all quantum states, represented as a list of n-iary \
(e.g. binary) numbers  *)
allstates[nlevels_, ntransmons_] := 
  Table[ IntegerDigits[ i, nlevels, ntransmons ] , { i , 0 , 
    nlevels^ntransmons - 1 } ];

(* Get all states with a total number of particles equal to \
nparticles *)
particlestates[ nlevels_, ntransmons_, nparticles_ ] :=  
  Select[ allstates[nlevels, ntransmons],  Total[#] == nparticles & ] ;
  
(* Give all states with nparticles, but as a more compact string *)
particlestatestr[ nlevels_, ntransmons_, nparticles_ ] :=  
  StringJoin/@ Map[ ToString, particlestates[ nlevels, ntransmons, nparticles ] , {2} ]
   
(* Get "allstates" but as a string, for pretty output *)
allstatestr[ nlevels_, ntransmons_ ] := 
  Table[StringJoin @@  
    ToString /@ IntegerDigits[ i, nlevels, ntransmons ] , { i , 0 , 
    nlevels^ntransmons - 1 } ];
    
(* Get all indices corresponding to states with nparticles *)
particleindices[ nlevels_, ntransmons_, nparticles_ ] :=  
  FromDigits[ #, nlevels] + 1  & /@ 
   particlestates[ nlevels, ntransmons, nparticles]  ;


(* Let operator "op" work on site "site", where we have nlevels=2 and \
nsites = #qubits. *)
qop[op_, site_, nsites_] := Module[{order},
   order = Log[2, Length[op]]; 
   KroneckerProduct[IdentityMatrix[2^(site - 1)], op, 
    IdentityMatrix[2^(nsites - site - (order - 1))]]
   ];
   
(* Create the unit vector e_a of dimension n. Counting goes from \
0 to (n-1) *)
Vec[a_, n_] := UnitVector[ n, a+1 ];

(* Creates a matrix all elements zero, except for [i,j], which is 1. *)
matrixBasis[i_, j_, n_] := Module[{tmp},
   tmp = ConstantArray[0, {n, n}];
   tmp[[i, j]] = 1;
   Return[tmp];
   ];

(* Solve time evolution starting in any state, for any Hamiltonian \
that is a function of t. Returns an InterpolatingFunction.  *)

evolveState[ ham_, startstate_, {t_, tmin_, tmax_} ] := Module[ {sol}, 
  sol = NDSolve[ {I  v'[t] == ham . v[t], v[0] ==  startstate}, 
    v, {t, tmin, tmax} ] ;
  v  /. sol // Flatten // First
  ];
  
  
(* New and improved version: Now outputs a *list* of InterpolatingFunctions, making plotting easier *)
evolveState2[ ham_, startstate_, {t_, tmin_, tmax_} ] := Module[ {vlist, rhs, lhs}, 
	vlist = Through[ Array[ v, Length[startstate] ][t] ];
	lhs = I D[ vlist , t ];
	rhs = ham.vlist;

	NDSolveValue[ And[
		LogicalExpand[lhs == rhs ], 
		LogicalExpand[(vlist /. t-> tmin ) == startstate] ],
		vlist,{t, tmin, tmax} ] 
  ];
  
minslices[ tmax_, w_, minsliceperosc_:8 ] := minsliceperosc * tmax * w / (2 Pi) ;


evolveExp[ ham_, {t_,t0_,tf_}, ntslices_, toffset_:0, onlyFinal_:False  ]:= Module[{us, dt, hamslice },

	us = ConstantArray[None, ntslices+1];
	us[[1]] = IdentityMatrix[Length[ham]];
	dt = (tf-t0)/ntslices;
	Do[
		hamslice = N[ham /. t-> ( i-1 + toffset)*dt];
		us[[i+1]] = MatrixExp[ -I hamslice dt  ].us[[i]];
	,{i,1,ntslices}];
	If[onlyFinal, Return[us[[-1]]], Return[us]];
];

(* Philosophy: Apply evolveExp over one period, then use MatrixPower to get the final unitary 
	Note: Period is 2\[Pi]/\[Omega] for angular velocity \[Omega]. *)
evolveExpFloquet[ ham_, t_,period_, nperiods_, ntslices_, toffset_:0, onlyFinal_:True  ]:= Module[{oneperiod},
	oneperiod = evolveExp[ ham, {t,0,period}, ntslices, toffset, True ];
	If[ onlyFinal, Return[MatrixPower[ oneperiod, nperiods ]], Table[ MatrixPower[ oneperiod, i ], {i,0,nperiods} ]];
];

(* Solve the time evolution of unitary propagator up to time tau  *)
evolveUnitary[ ham_, {t_, tmin_, tmax_} ] := 
  NDSolveValue[ { U'[t] == -I ham . U[t], 
    U[0] == IdentityMatrix[Length[ham]] }, U, {t, tmin, tmax} ];


(* new and improved version: Outputs a LIST of interpolatingfunctions. 
	Solve the time evolution of unitary propagator up to time tau  *)
evolveUnitary2[ ham_, {t_, tmin_, tmax_} ] := Module[{d, vlist, lhs, rhs},

	d = Length[ham];
	vlist = Through /@ Through[ Array[ v, {d,d} ][t]  ] ;
	lhs = I D[ vlist , t ];
	rhs = ham.vlist;

	NDSolveValue[ And[
		LogicalExpand[lhs == rhs ], 
		LogicalExpand[(vlist /. t-> tmin ) == IdentityMatrix[d] ] ],
		vlist,{t, tmin, tmax} ] 
];

(* Split vector of interpolationfunctions *)
interpolationElement[ ifun_, index_ ] := Module[ {pts},
  pts = Transpose[
    Append[ifun["Coordinates"], ifun["ValuesOnGrid"][[All, index]] ]] ;
  Interpolation[pts]
  ]
  
(* Make a tensor with dimensions (nfunctions, ntimesteps, 2) (last \
dimension is (X,Y) when plotting) from interpolatingfunctions. This \
is the format accepted by ListPlot. *)
makePlotGrid[ sol_ , skip_: 1 ] := Module[ {gridValues, gridCoords} ,
   gridValues = sol["ValuesOnGrid"];
   gridCoords = sol["Coordinates"];
    Table[ {gridCoords[[1, 1 ;; -1 ;; skip]], 
      gridValues[[1 ;; -1 ;; skip, i]] }\[Transpose], {i, 1, 
     Length[gridValues[[1]]]}]
   ];
   
(*******************  Newly added : ************************)

(* Gives indices of the p-particle sector (for n qubits), assuming \
particle-ordered eigenbasis (eg 0001, 0010, 0100, etc )  *)
pidx[p_, n_] := Module[{start, end},
   start = Sum[ Binomial[n, i], {i, 0, p - 1} ] + 1;
   end = start + Binomial[ n, p ] - 1;
   Range[2^n][[start ;; end]]
   ];



End[]; (* End the private context *)
EndPackage[];
