(* ::Package:: *)

(* Credits to:

Stackoverflow user Everett You. Source: \
http://mathematica.stackexchange.com/questions/54545/is-it-possible-\
to-define-a-new-plottheme

Chris1992 on https://mathematica.stackexchange.com/questions/5132/export-plots-to-latex?noredirect=1&lq=1

*)


BeginPackage[ "KoensPlotStyle`"]


Begin["`Private`"]

Themes`AddThemeRules["Koen",
	DefaultPlotStyle->Thread@Directive[ColorData[112,"ColorList"], Thick ],
	Background->None,
	Axes-> False,
	Frame->True,
	AxesStyle->Directive[AbsoluteThickness[1],Black,FontSize->14],FrameStyle->Directive[AbsoluteThickness[1],Black,FontSize->14],
	TicksStyle->Directive[Black,FontSize->12],
	FrameTicksStyle->Directive[Black,FontSize->12],
	GridLinesStyle->Directive[AbsoluteThickness[0.5],Opacity[0.5]],
	ImageSizeRaw->{{400},{400}},
	LabelStyle->Directive[Black,FontSize->12],
	PlotLegends-> Automatic
];


$PlotTheme="Koen";

End[]  (* End the private context *)

EndPackage[]

