(* ::Package:: *)

BeginPackage[ "KrawtchoukFunctions`"];
Needs["QuantumInformation`"];

(* Public section --- Everything defined here is directly evaluated. Any parameters referenced will turn black *)

(* OLD QUANTUM INFORMATION PART - NOW IMPORTED SEPERATELY ---->

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
allstatestr::usage = "allstatestr[nlevels, ntransmons] gives a string representing all fock states";
particleindices::usage = "particleindices[nlevels, ntransmons, nparticles] gives the indices of nparticle states";

qop; Vec; matrixBasis; evolveState; evolveUnitary; interpolationElement; makePlotGrid; pidx;

END OF OLD QINF PART *)


(* Define functions publicly *)

krawtchouk; hamk; kceigensystem; hamZ; ZcorrectionFromHam;  eigengateThreeExp; eigengateOneExp; eigengate;
kchouck; fi; qparticlestates; xsToIdx; ampl; kcstate; evecSlater; evecSlaterBinOrder; 



(* Private section. Any defined functions/modules should placed here *)
Begin["`Private`"];


(* OLD QINF FUNCTIONS BEGIN ----> 

(*  Get a list of all quantum states, represented as a list of n-iary \
(e.g. binary) numbers  *)

allstates[nlevels_, ntransmons_] := 
  Table[ IntegerDigits[ i, nlevels, ntransmons ] , { i , 0 , 
    nlevels^ntransmons - 1 } ];
(* Get all states with a total number of particles equal to \
nparticles *)

particlestates[ nlevels_, ntransmons_, nparticles_ ] :=  
  Select[ allstates[nlevels, ntransmons], 
   Total[#] == nparticles & ] ;
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
(* Create the unit vector e_a with total length n. Counting goes from \
0 to (n-1) *)

Vec[a_, n_] := Table[KroneckerDelta[a + 1, i], {i, 1, n}];
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
  ]

(* Solve the time evolution of unitary propagator up to time tau  *)

evolveUnitary[ ham_, {t_, tmin_, tmax_} ] := 
  NDSolveValue[ { U'[t] == -I ham . U[t], 
    U[0] == IdentityMatrix[Length[ham]] }, U, {t, tmin, tmax} ];

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
   
(* Gives indices of the p-particle sector (for n qubits), assuming \
particle-ordered eigenbasis  *)
pidx[p_, n_] := Module[{start, end},
   start = Sum[ Binomial[n, i], {i, 0, p - 1} ] + 1;
   end = start + Binomial[ n, p ] - 1;
   Range[2^n][[start ;; end]]
   ];
   
  End of old QINF functions *)



(************************* KRAWTCHOUK STUFF ******************************)

(* Definitions of krawtchouk chain *)

krawtchouk[nq_] := 1/2 Sqrt[  # (nq - #) ] & /@ Range[nq - 1] ;
hamk[nq_] := 
  Sum[ krawtchouk[nq][[i]] * qop[XXYY/2, i, nq] , {i, 1, nq - 1} ];

(* These are the eigenstates, as found by diagonalizing hamk[] *)

kceigensystem[ nq_ ] := Module[ {eval, evec, order, evecn, hameb} ,
   {eval, evec} = 
    Eigensystem[ 
     N[ hamk[nq]  + Sum[-300 qop[Z, i, nq] , {i, 1, nq} ] ] ];
   order = Ordering[eval];
   evec = evec[[order]];
   evecn = Chop[ Normalize /@ evec , 10^-10];
   
   hameb = evecn . hamk[nq] . evecn\[Transpose] // Simplify;
   eval = Chop[ Diagonal[hameb] ];
   {eval, evecn}
   ];

hamZ[nq_] := 
  1/2 Sum[   ((nq - 1)/2 - i) qop[Z , i + 1, nq ], {i, 0, nq - 1}] ;
ZcorrectionFromHam[nq_] := MatrixExp[ -I N[ hamZ[nq] Pi/2 ] ];

eigengateThreeExp[nq_] :=  
  ZcorrectionFromHam[nq] . MatrixExp[ -I N[ hamk[nq] Pi/2] ] . 
   ZcorrectionFromHam[nq];
eigengateOneExp[nq_] := 
  MatrixExp[ -I N[ (hamk[nq] + hamZ[nq])/Sqrt[2]  Pi ] ];

(* This is the action of the eigengate by explicit MatrixExp[ hamk[] \
] *)
eigengate[nq_] := Module[ {Hexp, Zcorrection} ,
   Hexp = MatrixExp[ -I N[ hamk[nq] Pi/2] ];
   Zcorrection = 
    KroneckerProduct @@ 
     Table[ DiagonalMatrix[{1, (-I)^j}], {j, 0, nq - 1}];
   Chop[  Zcorrection.Hexp.Zcorrection ]
   ];



(*
For the next functions, the number of qubits is n+1.

Now, we want to construct the many-particle eigenstates, using the slater determinant.
Input: {k_1, k_2, ... k_q} for q eigenparticles present.

| k_1, ... k_q >  = Sum_{x1 < x2 < ... < xq }  ampl[  {x1, ... x_q}, {k_1, ... k_q} ]   | x1, ... x_q >
*)


(* Definition of Krawtchouk matrix *)
Clear[kchouck, n, c, w, fi];
kchouck[ n_, k_, x_ ] := 
  Sum[ (-1)^i  Binomial[ x, i ] Binomial[ n - x, k - i] , { i, 0 , 
    k } ];
c[ n_, k_, x_ ] := Sqrt[ Pochhammer[-n, k] / (-1)^k / k!  ];
w[ n_, k_, x_ ] := 1/2^n * Binomial[ n , x ];
fi[n_, k_, x_] := 
  Sqrt[  Binomial[ n , x ] / Binomial[n, k]  / 
     2^n ] Sum[ (-1)^i  Binomial[ x, i ] Binomial[ n - x, k - i] , { 
     i, 0 , k } ];

(* Enumerate all n-qubit states with q particles, such as { {0,1,4}, ... \
} *)
qparticlestates[ n_, q_ ] := 
  Select[ Subsets[Range[0, n]], Length[#] == q & ];


(* Given a "particle presence" vector |1,0,1,0,0...>, give the \
"binary counting" vector index (e.g. [[16]]) *)
xsToIdx[ xs_ ] := Sum[ 2^xs[[i]] , {i, 1, Length[xs] } ];

(* Amplitude of the state | xs > in the state | ks > *)

ampl[ n_, xs_, ks_ ] := Module[{q},
  If[ Length[xs] < 1, Return[ 1 ] ];
  q = Length[ xs ];
  Return[ 
   Det[ Table[ 
     fi[n, xs[[idx]], ks[[idk]] ], {idx, 1, q}, {idk, 1, q} ] ] ];
  ]

(* Build the state | ks > in the "binary counting" representation *)

kcstate[ n_, ks_  ] := Module[{relevantxs},
  If[ Length[ks] < 1 , Return[ Vec[0, 2^(n + 1) ] ] ];
  relevantxs = qparticlestates[ n, Length[ks ] ];
  Sum[   ampl[ n, relevantxs[[i]], ks ] Vec[ 
     xsToIdx[ relevantxs[[i]] ], 2^(n + 1) ] , {i, 1, 
    Length[ relevantxs ] } ]
  ]

(* Make the whole Krawtchouk basis. Ordering according to the \
spectrum (first particle number, then energy) *)

evecSlater[n_] := Module[ {allks},
   allks = Subsets[ Range[0, n] ];
   Table[ kcstate[ n, allks[[i]] ], {i, 1, Length[allks] } ]
   ];

(* Same as previous, but ordered according to binary ordering of \
inputs *)
evecSlaterBinOrder[n_] := Module[ {xs},
  Table[
   (* We turn the binary number i into the list of indices indicating \
the 1's *)
   
   xs = Position[ Reverse@IntegerDigits[ i, 2, n + 1 ] , 1 ] - 1 // 
     Flatten;
   kcstate[n, xs], {i, 0, 2^(n + 1) - 1} ]
  ]



End[]; (* End the private context *)

EndPackage[];
