(* ::Package:: *)

(* This is the package GRAPE.wl (import as Needs[GRAPE`] ) which contains my personal functions for GRAPE optimization. *)

 BeginPackage[ "GRAPEKoen`"];

(* Public section --- Everything defined here is directly evaluated. Any parameters referenced will turn black *)

(* Define all functions / modules *)
fidelity::usage = "placeholder";
fidNoSquare::usage = "placeholder";
forwardPass::usage = "placeholder";
makeUnitary::usage = "placeholder";
gradients::usage = "placeholder";
stupidgrad::usage = "placeholder";
doGrapeOptimization::usage = "placeholder";
FourierGrape::usage = "placeholder";





(* ======================================================================================================= *)
(* Private section. Any defined functions/modules should placed here *)
Begin["`Private`"];

(* Calculate the Fidelity (square of the Hilbert-Schmidt inner product) *)
fidelity[ ua_, ub_ ]:= Abs[ Tr[ConjugateTranspose[ua].ub]/Length[ua]  ]^2;
fidNoSquare[ua_,ub_]:= Abs[ Tr[ConjugateTranspose[ua].ub] ] / Length[ua];


(* 
Find the simulated unitary based on ALL timeslices 
forwardPass[ baseham, fields, u, dt, ntslices ]
*)
forwardPass[ baseham_, fields_, u_, dt_, ntslices_  ] := Module[{prod, slice},
	(* u has shape (nfields, nslices) *)
	prod  =  IdentityMatrix[ Length[baseham] ];
	For[ slice = 1, slice<= ntslices, slice ++,
		prod = MatrixExp[ -I ( baseham + Sum[ u[[f,slice]] * fields[[f]] , {f,1,Length[fields]}  ] ) dt ] . prod
	];
	Return[prod]; (* return value *)
];


(* Find the unitary corresponding to the given timeslice "slice" *)
makeUnitary[ baseham_, fields_, u_, dt_, slice_ ]  := MatrixExp[ -I ( baseham + 
	Sum[ u[[f,slice]] * fields[[f]] , {f,1,Length[fields]}  ] ) dt ];


(* 
Find gradient of the fidelity (squared) with respect to u[f,slice]. Returns a matrix in the same shape as u. 
gradients[  utarget, baseham, fields, u, dt, ntslices  ]
*)
gradients[  utarget_, baseham_, fields_, u_, dt_, ntslices_  ] := Module [ {unitariesPerSlice, du, Xj, Pjdag, slice, f, newUnitary},

	(* Pre-calculate all the (slow!) matrix exponentials *)
	unitariesPerSlice = Table[ makeUnitary[ baseham, fields, u, dt, slice ], {slice, 1,ntslices} ];


	(* For every iteration, set these initial values.... *)
	du = ConstantArray[ 0, {Length[fields], ntslices}  ];
	Xj  =  unitariesPerSlice[[1]]; (* take slice 1 *)
	Pjdag = ConjugateTranspose[  utarget ] .  Apply[ Dot,Reverse@ unitariesPerSlice[[2;;]] ] ;

	(* Loop over all timeslices *)
	For[ slice = 1, slice<= ntslices, slice ++,
		(* Loop over all fields *)
		For[ f = 1, f <= Length[fields], f++,

			du[[f, slice]] = -2 Re[   I dt Tr[ Pjdag . fields[[f]] .Xj ]  * (Conjugate@ Tr[ Pjdag.Xj ]) ] / Length[baseham]^2;

		];

		(* Prepare the next timeslice: *)
		If[ slice+1 <= ntslices ,
			newUnitary = makeUnitary[ baseham, fields, u, dt, slice+1 ];
			Xj = newUnitary.Xj;
			Pjdag = Pjdag.ConjugateTranspose[newUnitary];
		];

	];  (* end of timeslice loop *)

	Return[du];

] (* End of gradients function *)


(* Very inefficient way to calculate gradients... but works for large time-slices! *)
stupidgrad[utarget_, baseham_, fields_, u_, dt_, ntslices_  ] := Module[{findiffstep, chng, f2, f1},
	findiffstep = 0.01;
	Table[
		Table[
			chng = ConstantArray[0, Dimensions[u] ];
			chng[[fslice,tslice]] = findiffstep;

			f2 = fidelity[ forwardPass[ baseham, fields, u+chng, dt, ntslices  ], utarget ];
			f1 = fidelity[ forwardPass[ baseham, fields, u-chng, dt, ntslices  ], utarget ];
			(f2 - f1 )/(2findiffstep),
		{tslice, 1, ntslices}],
	{fslice, 1, Length[u]}] 
];


(*
Function that automatically optimizes u based on the input.
Verbose: set to 0 to see initialization settings. Set to -1 for no output.

doGrapeOptimization[  utarget, baseham, fields, ntslices, dt, learnstep, nruns, verbose, uinit] 
*)
doGrapeOptimization[  utarget_, baseham_, fields_, ntslices_, dt_, learnstep_, nruns_, verbose_: -1, uinit_: None] := Module[ {u, du,run, stdrange, relchange},
	(* 
	Initialize u, where uinit can be either 1) the complete matrix u or 2) an initialization range *per field* or 3)  an initialization range for all parameters. 
	*)
	stdrange = .3;
	Which[
(Length[uinit] == Length[fields] && Length[uinit[[1]] ] == ntslices ),
		If[verbose>= 0,  Print["Setting u = uinit"]]; 
		u = uinit;,
(Length[uinit] == Length[fields] ), 
		If[verbose>= 0,Print["Random u, separate range per field"]]; 
		u = uinit *RandomReal[ {-1, 1}, {Length[fields], ntslices} ] ;,
(NumberQ[uinit]), 
		If[verbose>= 0,Print["Generating initial u in range ", uinit ]];
		 u = RandomReal[ {-uinit, uinit}, {Length[fields], ntslices} ]; ,
(True), 
		 If[verbose>= 0,Print["Generating initial u in standard range ", stdrange]]; 
		 u = RandomReal[ {-stdrange, stdrange}, {Length[fields], ntslices} ]; 
	];

(* 
	Optimize that gate!
 *)
	If[verbose>= 1, Print["run \t error \t range of u \t \t \t \t \t range of du \t \t \t rel. change" ]];
	For[ run = 1, run <= nruns, run++,
		du = gradients[ utarget, baseham, fields, u, dt, ntslices ] ;
		u = u + learnstep * du;

		If[verbose >= 1 &&  Mod[run, verbose] == 0, 
			relchange = -learnstep * du / u;
			Print[ run, "\t",1-fidelity[ utarget, forwardPass[ baseham, fields, u, dt, ntslices ] ], "\t ", {Min[u], Max[u]}, "\t \t",  {Min[du],Max[du]}, "\t \t", {Min[relchange], Max[relchange]}  ];
		];
	]; (* End of For loop over nruns *)

	Return[u]
];

(* According to documentation, the normal Fourier does the following: 1/Sqrt[n]\!\(\*
UnderoverscriptBox["\[Sum]", 
RowBox[{
StyleBox["r", "TI"], "=", 
StyleBox["1", "TR"]}], 
StyleBox["n", "TI"],
LimitsPositioning->True]\(\*
SubscriptBox[
StyleBox["u", "TI"], 
StyleBox["r", "TI"]]\*
SuperscriptBox[
StyleBox["e", "TI"], 
RowBox[{
StyleBox["2", "TR"], 
StyleBox["\[Pi]", "TR"], 
StyleBox[" ", "TR"], 
StyleBox["i", "TI"], 
RowBox[{"(", 
RowBox[{
StyleBox["r", "TI"], "-", 
StyleBox["1", "TR"]}], ")"}], 
RowBox[{
RowBox[{"(", 
RowBox[{
StyleBox["s", "TI"], "-", 
StyleBox["1", "TR"]}], ")"}], "/", 
StyleBox["n", "TI"]}]}]] . \ \ The\ values\ of\ r\ represent\ time\)\), and we pick r-1/2 instead of r-1. The frequency domain (s) remains unchanged. *)

FourierGrape[ u_ ] := Module[{n},
n = Length[u];
1/Sqrt[n] * Table[    
Sum[ u[[r]]*  Exp[ 2 Pi I (r-1/2)(s-1) / n ], {r,1,n}]
,{s,1,n}]
];










End[]; (* End the private context *)

EndPackage[];

