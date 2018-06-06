(* Wolfram Language Package *)

(* Created by the Wolfram Workbench 10-Apr-2015 *)

BeginPackage["BarrierOption`"]
(* Exported symbols added here with SymbolName::usage *) 

MakeGoodCorrelationMatrix::usage="MakeGoodCorrelationMatrix[badCorrMat,flag] returns an acceptable correlation matrix (positive semi-definite)
	which is close to the input one. Close means that the sum of the squares of the difference of all the elements is minimized.
	The input matrix is assumed to be symmetric.

	flag = 1 uses the hypersphere decomposition
	flag != 1 uses the spectral decomposition as detailed in the paper by R. Rebonato and P. Jackel."
bp::usage="Helper function"
bpp::usage="Helper function"
outCoeff::usage

(*defgGeneric::usage=""*)
phi::usage="Utility function"
smallPhi::usage=""
bigPhi::usage=""
localDot::usage=""

barrierOption::usage="Returns the value of a barrier option"
doubleKnockOutOutsideBarrierOptionValue::usage="Returns the value of a double knock-out outside barrier option."
doubleKnockInOutsideBarrierOptionValue::usage="Returns the value of a double knock-in outside barrier option."
singleUpAndOutOutsideBarrierOptionValue::usage"Returns the value of a single up-and-out outside barrier option."
singleUpAndInOutsideBarrierOptionValue::usage"Returns the value of a single up-and-in outside barrier option."
singleDownAndOutOutsideBarrierOptionValue::usage"Returns the value of a single down-and-out outside barrier option."
singleDownAndInOutsideBarrierOptionValue::usage"Returns the value of a single down-and-in outside barrier option."

Begin["`Private`"]
(* Implementation of the package *)

MakeGoodCorrelationMatrix[badCorrMat_, flag_] :=
    Which[
      flag == 1,
      Module[{
				LocalN = Length[badCorrMat], LocalTheta, LocalB, LocalOutput, LocalSolution},
        LocalTheta = Outer[Unique[] &, Range[LocalN], Range[LocalN - 1], 1, 1];
        LocalB = Outer[If[#2 < LocalN, Cos[LocalTheta[[#1, #2]]] (Times @@ Sin[LocalTheta[[#1, 1 ;; #2 - 1]]]),
          Times @@ Sin[LocalTheta[[#1, 1 ;; #2 - 1]]]] &, Range[LocalN], Range[LocalN], 1, 1] // Chop;
        LocalOutput = Dot[LocalB, Transpose[LocalB]] // Chop;
        LocalSolution = Minimize[Total[Flatten[(LocalOutput - badCorrMat)^2]], Flatten[LocalTheta]] // Chop;
        LocalOutput /. LocalSolution[[2]]

      ],
      True,
      Module[{
				LocalLambdap = (Max[10^6 $MachineEpsilon, #] & /@ Eigenvalues[badCorrMat]),
				LocalS = Eigenvectors[badCorrMat] // Transpose,
				LocalT, LocalB, LocalOutput},
				LocalT = DiagonalMatrix[Dot[#^2, LocalLambdap]^(-1) & /@ LocalS] ;
        LocalB = Dot[Sqrt[LocalT], LocalS, Sqrt[DiagonalMatrix[LocalLambdap]]];

        LocalOutput = Dot[LocalB, Transpose[LocalB]] ;
        LocalOutput = 1 / 2 (LocalOutput + Transpose[LocalOutput]) (* just to make sure the output is symmetric *)

      ]
    ]


daysYear = 365.242;

bp[low_, high_, n_] = 2 (high - low) n;
bpp[low_, high_, n_] = 2 high -  bp[low, high, n];

defgGeneric[nontterm_, tterm_, t_, lastterm_, den_, additionalterm_] = (nontterm + tterm t + lastterm) / den + additionalterm;

(*localDot[x_, y_] /; And[VectorQ[x], VectorQ[y]] := Total[x y];
localDot[x_, y_] /; And[MatrixQ[x], VectorQ[y]] := Map[Total[# y] &, x];
localDot[x_, y_] /; And[VectorQ[x], MatrixQ[y]] := Map[Total[x #] &, y];
localDot[x_, y_] /; And[MatrixQ[x], MatrixQ[y]] := Outer[Total[#1 #2] &, x, Transpose[y], 1, 1]
SetAttributes[localDot, NumericFunction];*)

pdfMultiNormalDistribution[corrMat_, x_] = 1 / ((2 Pi)^(Length[corrMat] / 2) * Sqrt[Det[corrMat]]) Exp[-1/2 Dot[x, Dot[Inverse[corrMat], x]]];
SetAttributes[pdfMultiNormalDistribution, Union[Attributes[pdfMultiNormalDistribution], {NumericFunction}]];
smallPhi[mean_, covMat_, x_] := 
	Module[{vols, corrMat}, 
		vols = Sqrt[Diagonal[covMat]];
		corrMat = covMat / Outer[#1 #2 &, vols ,vols, 1, 1]; 
		pdfMultiNormalDistribution[corrMat, (x - mean) / vols]
	];
SetAttributes[smallPhi, Union[Attributes[smallPhi], {NumericFunction}]];
	
(*phi = Function[{mean, cov}, vols=Sqrt[Diagonal[cov]]; CDF[MultinormalDistribution[ConstantArray[0, Length[cov]], cov / Outer[#1 #2 &, vols ,vols, 1, 1]], (# - mean) / vols] &];*)
phi[mean_, covMat_, x_] := 
	Module[{len, vars, localPhi}, 
		len = Length[x];
		vars = Unique["z"] & /@ Range[len];
		localPhi = smallPhi[mean, covMat, vars];
		NIntegrate@@{localPhi, Sequence @@ Transpose[{vars, ConstantArray[-Infinity, len], x}]}
	]
SetAttributes[phi, NumericFunction];


outCoeff[cp_, {sUnd_, sYield_, sVol_}, {bUnd_, bYield_, bVol_}, corr_, intRate_, {t1_, t2_, bLow_, bHigh_}, strike_, t_, ts_, n_] :=
    Module[ {
            a, bH, bL, mu1, mu2, 
            gMat, d1p, d2p, d1pp, d2pp, e1p, e2p, e1pp, e2pp, f1p, f2p, f1pp, f2pp, g1, g2, h1, h2,
            wsn, wkn},
        a = Log[strike / sUnd];
        bL = Log[bLow / bUnd];
        bH = Log[bHigh / bUnd];
        mu1 = intRate - sYield - 1/2 sVol^2;
        mu2 = intRate - bYield - 1/2 bVol^2;
        gMat = {{1, -cp corr Sqrt[t2/t], -cp corr Sqrt[t1/t]}, {-cp corr Sqrt[t2/t], 1, Sqrt[t1/t2]}, {-cp corr Sqrt[t1/t], Sqrt[t1/t2], 1}};
        d1p = defgGeneric[-a, mu1, t, corr sVol/bVol bp[bL, bH, n], sVol Sqrt[t], sVol Sqrt[t]];
        d2p = defgGeneric[-a, mu1, t, corr sVol/bVol bp[bL, bH, n], sVol Sqrt[t], 0];
        d1pp = defgGeneric[-a, mu1, t, corr sVol/bVol bpp[bL, bH, n], sVol Sqrt[t], sVol Sqrt[t]];
        d2pp = defgGeneric[-a, mu1, t, corr sVol/bVol bpp[bL, bH, n], sVol Sqrt[t], 0];
        e1p = defgGeneric[bH, -mu2, t2, -bp[bL, bH, n], bVol Sqrt[t2], -corr sVol Sqrt[t2]];
        e2p = defgGeneric[bH, -mu2, t2, -bp[bL, bH, n], bVol Sqrt[t2], 0];
        e1pp = defgGeneric[bH, -mu2, t2, -bpp[bL, bH, n], bVol Sqrt[t2], -corr sVol Sqrt[t2]];
        e2pp = defgGeneric[bH, -mu2, t2, -bpp[bL, bH, n], bVol Sqrt[t2], 0];
        f1p = defgGeneric[bL, -mu2, t2, -bp[bL, bH, n], bVol Sqrt[t2], -corr sVol Sqrt[t2]];
        f2p = defgGeneric[bL, -mu2, t2, -bp[bL, bH, n], bVol Sqrt[t2], 0];
        f1pp = defgGeneric[bL, -mu2, t2, -bpp[bL, bH, n], bVol Sqrt[t2], -corr sVol Sqrt[t2]];
        f2pp = defgGeneric[bL, -mu2, t2, -bpp[bL, bH, n], bVol Sqrt[t2], 0];
        g1 = defgGeneric[bH, -mu2, t1, 0, bVol Sqrt[t1], -corr sVol Sqrt[t1]];
        g2 = defgGeneric[bH, -mu2, t1, 0, bVol Sqrt[t1], 0];
        h1 = defgGeneric[bL, -mu2, t1, 0, bVol Sqrt[t1], -corr sVol Sqrt[t1]];
        h2 = defgGeneric[bL, -mu2, t1, 0, bVol Sqrt[t1], 0];
        wsn = Exp[(mu2 + corr sVol bVol) bp[bL, bH, n] / bVol^2] * 
                (CDF[MultinormalDistribution[{0, 0, 0}, gMat], {cp d1p, e1p, g1}] -
                 CDF[MultinormalDistribution[{0, 0, 0}, gMat], {cp d1p, f1p, g1}] -
                 CDF[MultinormalDistribution[{0, 0, 0}, gMat], {cp d1p, e1p, h1}] +
                 CDF[MultinormalDistribution[{0, 0, 0}, gMat], {cp d1p, f1p, h1}]) -
              Exp[(mu2 + corr sVol bVol) bpp[bL, bH, n] / bVol^2] * 
                (CDF[MultinormalDistribution[{0, 0, 0}, gMat], {cp d1pp, e1pp, g1}] -
                 CDF[MultinormalDistribution[{0, 0, 0}, gMat], {cp d1pp, f1pp, g1}] -
                 CDF[MultinormalDistribution[{0, 0, 0}, gMat], {cp d1pp, e1pp, h1}] +
                 CDF[MultinormalDistribution[{0, 0, 0}, gMat], {cp d1pp, f1pp, h1}]);
        wkn = Exp[mu2  bp[bL, bH, n] / bVol^2] * 
                (CDF[MultinormalDistribution[{0, 0, 0}, gMat], {cp d2p, e2p, g2}] -
                 CDF[MultinormalDistribution[{0, 0, 0}, gMat], {cp d2p, f2p, g2}] -
                 CDF[MultinormalDistribution[{0, 0, 0}, gMat], {cp d2p, e2p, h2}] +
                 CDF[MultinormalDistribution[{0, 0, 0}, gMat], {cp d2p, f2p, h2}]) -
              Exp[mu2 bpp[bL, bH, n] / bVol^2] * 
                (CDF[MultinormalDistribution[{0, 0, 0}, gMat], {cp d2pp, e2pp, g2}] -
                 CDF[MultinormalDistribution[{0, 0, 0}, gMat], {cp d2pp, f2pp, g2}] -
                 CDF[MultinormalDistribution[{0, 0, 0}, gMat], {cp d2pp, e2pp, h2}] +
                 CDF[MultinormalDistribution[{0, 0, 0}, gMat], {cp d2pp, f2pp, h2}]);
        {gMat, wsn, wkn}
    ]


barrierOption[cp_, {sUnd_, sYield_, sVol_}, {bUnd_, bYield_, bVol_}, corr_, intRate_, {t1_, t2_, bLow_, bHigh_}, strike_, t_, ts_, {nMin_, nMax_}] :=
    Module[ {
            a, bH, bL, mu1, mu2, 
            gMat, d1p, d2p, d1pp, d2pp, e1p, e2p, e1pp, e2pp, f1p, f2p, f1pp, f2pp, g1, g2, h1, h2,
            wsn, wkn, ws, wk, 
            value
        },
        a = Log[strike / sUnd];
        bL = Log[bLow / bUnd];
        bH = Log[bHigh / bUnd];
        mu1 = intRate - sYield - 1/2 sVol^2;
        mu2 = intRate - bYield - 1/2 bVol^2;
        gMat = {{1, -cp corr Sqrt[t2/t], -cp corr Sqrt[t1/t]}, {-cp corr Sqrt[t2/t], 1, Sqrt[t1/t2]}, {-cp corr Sqrt[t1/t], Sqrt[t1/t2], 1}};
        d1p = defgGeneric[-a, mu1, t, corr sVol/bVol bp[bL, bH, n], sVol Sqrt[t], sVol Sqrt[t]];
        d2p = defgGeneric[-a, mu1, t, corr sVol/bVol bp[bL, bH, n], sVol Sqrt[t], 0];
        d1pp = defgGeneric[-a, mu1, t, corr sVol/bVol bpp[bL, bH, n], sVol Sqrt[t], sVol Sqrt[t]];
        d2pp = defgGeneric[-a, mu1, t, corr sVol/bVol bpp[bL, bH, n], sVol Sqrt[t], 0];
        e1p = defgGeneric[bH, -mu2, t2, -bp[bL, bH, n], bVol Sqrt[t2], -corr sVol Sqrt[t2]];
        e2p = defgGeneric[bH, -mu2, t2, -bp[bL, bH, n], bVol Sqrt[t2], 0];
        e1pp = defgGeneric[bH, -mu2, t2, -bpp[bL, bH, n], bVol Sqrt[t2], -corr sVol Sqrt[t2]];
        e2pp = defgGeneric[bH, -mu2, t2, -bpp[bL, bH, n], bVol Sqrt[t2], 0];
        f1p = defgGeneric[bL, -mu2, t2, -bp[bL, bH, n], bVol Sqrt[t2], -corr sVol Sqrt[t2]];
        f2p = defgGeneric[bL, -mu2, t2, -bp[bL, bH, n], bVol Sqrt[t2], 0];
        f1pp = defgGeneric[bL, -mu2, t2, -bpp[bL, bH, n], bVol Sqrt[t2], -corr sVol Sqrt[t2]];
        f2pp = defgGeneric[bL, -mu2, t2, -bpp[bL, bH, n], bVol Sqrt[t2], 0];
        g1 = defgGeneric[bH, -mu2, t1, 0, bVol Sqrt[t1Local], -corr sVol Sqrt[t1]];
        g2 = defgGeneric[bH, -mu2, t1, 0, bVol Sqrt[t1Local], 0];
        h1 = defgGeneric[bL, -mu2, t1, 0, bVol Sqrt[t1Local], -corr sVol Sqrt[t1]];
        h2 = defgGeneric[bL, -mu2, t1, 0, bVol Sqrt[t1Local], 0];
        wsn = Exp[(mu2 + corr sVol bVol) bp[bL, bH, n] / bVol^2] * 
                (CDF[MultinormalDistribution[{0, 0, 0}, gMat], {cp d1p, e1p, g1}] -
                 CDF[MultinormalDistribution[{0, 0, 0}, gMat], {cp d1p, f1p, g1}] -
                 CDF[MultinormalDistribution[{0, 0, 0}, gMat], {cp d1p, e1p, h1}] +
                 CDF[MultinormalDistribution[{0, 0, 0}, gMat], {cp d1p, f1p, h1}]) -
              Exp[(mu2 + corr sVol bVol) bpp[bL, bH, n] / bVol^2] * 
                (CDF[MultinormalDistribution[{0, 0, 0}, gMat], {cp d1pp, e1pp, g1}] -
                 CDF[MultinormalDistribution[{0, 0, 0}, gMat], {cp d1pp, f1pp, g1}] -
                 CDF[MultinormalDistribution[{0, 0, 0}, gMat], {cp d1pp, e1pp, h1}] +
                 CDF[MultinormalDistribution[{0, 0, 0}, gMat], {cp d1pp, f1pp, h1}]);
        wkn = Exp[mu2  bp[bL, bH, n] / bVol^2] * 
                (CDF[MultinormalDistribution[{0, 0, 0}, gMat], {cp d2p, e2p, g2}] -
                 CDF[MultinormalDistribution[{0, 0, 0}, gMat], {cp d2p, f2p, g2}] -
                 CDF[MultinormalDistribution[{0, 0, 0}, gMat], {cp d2p, e2p, h2}] +
                 CDF[MultinormalDistribution[{0, 0, 0}, gMat], {cp d2p, f2p, h2}]) -
              Exp[mu2 bpp[bL, bH, n] / bVol^2] * 
                (CDF[MultinormalDistribution[{0, 0, 0}, gMat], {cp d2pp, e2pp, g2}] -
                 CDF[MultinormalDistribution[{0, 0, 0}, gMat], {cp d2pp, f2pp, g2}] -
                 CDF[MultinormalDistribution[{0, 0, 0}, gMat], {cp d2pp, e2pp, h2}] +
                 CDF[MultinormalDistribution[{0, 0, 0}, gMat], {cp d2pp, f2pp, h2}]);
        ws = Sum[Evaluate[wsn], {n, nMin, nMax}];
        wk = Sum[Evaluate[wkn], {n, nMin, nMax}];
        value = Exp[-intRate (ts - t)] cp (sUnd Exp[-sYield t] ws - strike Exp[-intRate t] wk)
    ]


(* Quantities used in limiting cases *)
$zeroTime = 10^(-3);
$zeroBarrier = 10^(-4);
$infiniteBarrier = 10^4;
$unitCorrelation = 0.999;

doubleKnockOutOutsideBarrierOptionValue[cp_, {sUnd_, sYield_, sVol_}, {bUnd_, bYield_, bVol_}, corr_, intRate_, {t1_, t2_, bLow_, bHigh_}, strike_, t_, ts_, {nMin_, nMax_}] /; Or[bUnd < bLow, bUnd > bHigh] :=
    0
doubleKnockOutOutsideBarrierOptionValue[cp_, {sUnd_, sYield_, sVol_}, {bUnd_, bYield_, bVol_}, corr_, intRate_, {t1_, t2_, bLow_, bHigh_}, strike_, t_, ts_, {nMin_, nMax_}] :=
    Module[{t1Local, n, 
            a, bH, bL, mu1, mu2, 
            gMat, d1p, d2p, d1pp, d2pp, e1p, e2p, e1pp, e2pp, f1p, f2p, f1pp, f2pp, g1, g2, h1, h2,
            wsn, wkn, ws, wk, 
            value
        },
        t1Local = Max[t1, 10^(-3)(*QuantityMagnitude[UnitConvert[Quantity[1, "Hours"], "Years"]]*)]; (* Takes care of t1=0 in g1, g2, h1, h2 *)
        
        a = Log[strike / sUnd];
        bL = Log[bLow / bUnd];
        bH = Log[bHigh / bUnd];
        mu1 = intRate - sYield - 1/2 sVol^2;
        mu2 = intRate - bYield - 1/2 bVol^2;
        gMat = {{1, -cp corr Sqrt[t2/t], -cp corr Sqrt[t1/t]}, {-cp corr Sqrt[t2/t], 1, Sqrt[t1/t2]}, {-cp corr Sqrt[t1/t], Sqrt[t1/t2], 1}};
        gMat = If[Not[PositiveDefiniteMatrixQ[gMat]], MakeGoodCorrelationMatrix[gMat, 1], gMat];
        d1p = defgGeneric[-a, mu1, t, corr sVol/bVol bp[bL, bH, n], sVol Sqrt[t], sVol Sqrt[t]];
        d2p = defgGeneric[-a, mu1, t, corr sVol/bVol bp[bL, bH, n], sVol Sqrt[t], 0];
        d1pp = defgGeneric[-a, mu1, t, corr sVol/bVol bpp[bL, bH, n], sVol Sqrt[t], sVol Sqrt[t]];
        d2pp = defgGeneric[-a, mu1, t, corr sVol/bVol bpp[bL, bH, n], sVol Sqrt[t], 0];
        e1p = defgGeneric[bH, -mu2, t2, -bp[bL, bH, n], bVol Sqrt[t2], -corr sVol Sqrt[t2]];
        e2p = defgGeneric[bH, -mu2, t2, -bp[bL, bH, n], bVol Sqrt[t2], 0];
        e1pp = defgGeneric[bH, -mu2, t2, -bpp[bL, bH, n], bVol Sqrt[t2], -corr sVol Sqrt[t2]];
        e2pp = defgGeneric[bH, -mu2, t2, -bpp[bL, bH, n], bVol Sqrt[t2], 0];
        f1p = defgGeneric[bL, -mu2, t2, -bp[bL, bH, n], bVol Sqrt[t2], -corr sVol Sqrt[t2]];
        f2p = defgGeneric[bL, -mu2, t2, -bp[bL, bH, n], bVol Sqrt[t2], 0];
        f1pp = defgGeneric[bL, -mu2, t2, -bpp[bL, bH, n], bVol Sqrt[t2], -corr sVol Sqrt[t2]];
        f2pp = defgGeneric[bL, -mu2, t2, -bpp[bL, bH, n], bVol Sqrt[t2], 0];
        g1 = defgGeneric[bH, -mu2, t1Local, 0, bVol Sqrt[t1Local], -corr sVol Sqrt[t1Local]];
        g2 = defgGeneric[bH, -mu2, t1Local, 0, bVol Sqrt[t1Local], 0];
        h1 = defgGeneric[bL, -mu2, t1Local, 0, bVol Sqrt[t1Local], -corr sVol Sqrt[t1Local]];
        h2 = defgGeneric[bL, -mu2, t1Local, 0, bVol Sqrt[t1Local], 0];
        wsn = Exp[(mu2 + corr sVol bVol) bp[bL, bH, n] / bVol^2] * 
                (CDF[MultinormalDistribution[{0, 0, 0}, gMat], {cp d1p, e1p, g1}] -
                 CDF[MultinormalDistribution[{0, 0, 0}, gMat], {cp d1p, f1p, g1}] -
                 CDF[MultinormalDistribution[{0, 0, 0}, gMat], {cp d1p, e1p, h1}] +
                 CDF[MultinormalDistribution[{0, 0, 0}, gMat], {cp d1p, f1p, h1}]) -
              Exp[(mu2 + corr sVol bVol) bpp[bL, bH, n] / bVol^2] * 
                (CDF[MultinormalDistribution[{0, 0, 0}, gMat], {cp d1pp, e1pp, g1}] -
                 CDF[MultinormalDistribution[{0, 0, 0}, gMat], {cp d1pp, f1pp, g1}] -
                 CDF[MultinormalDistribution[{0, 0, 0}, gMat], {cp d1pp, e1pp, h1}] +
                 CDF[MultinormalDistribution[{0, 0, 0}, gMat], {cp d1pp, f1pp, h1}]);
        wkn = Exp[mu2  bp[bL, bH, n] / bVol^2] * 
                (CDF[MultinormalDistribution[{0, 0, 0}, gMat], {cp d2p, e2p, g2}] -
                 CDF[MultinormalDistribution[{0, 0, 0}, gMat], {cp d2p, f2p, g2}] -
                 CDF[MultinormalDistribution[{0, 0, 0}, gMat], {cp d2p, e2p, h2}] +
                 CDF[MultinormalDistribution[{0, 0, 0}, gMat], {cp d2p, f2p, h2}]) -
              Exp[mu2 bpp[bL, bH, n] / bVol^2] * 
                (CDF[MultinormalDistribution[{0, 0, 0}, gMat], {cp d2pp, e2pp, g2}] -
                 CDF[MultinormalDistribution[{0, 0, 0}, gMat], {cp d2pp, f2pp, g2}] -
                 CDF[MultinormalDistribution[{0, 0, 0}, gMat], {cp d2pp, e2pp, h2}] +
                 CDF[MultinormalDistribution[{0, 0, 0}, gMat], {cp d2pp, f2pp, h2}]);
        ws = Sum[Evaluate[wsn], {n, nMin, nMax}];
        wk = Sum[Evaluate[wkn], {n, nMin, nMax}];
        value = Exp[-intRate (ts - t)] cp (sUnd Exp[-sYield t] ws - strike Exp[-intRate t] wk)
    ]


doubleKnockInOutsideBarrierOptionValue[cp_, {sUnd_, sYield_, sVol_}, {bUnd_, bYield_, bVol_}, corr_, intRate_, {t1_, t2_, bLow_, bHigh_}, strike_, t_, ts_, {nMin_, nMax_}] /; Or[bUnd < bLow, bUnd > bHigh] :=
    0
doubleKnockInOutsideBarrierOptionValue[cp_, {sUnd_, sYield_, sVol_}, {bUnd_, bYield_, bVol_}, corr_, intRate_, {t1_, t2_, bLow_, bHigh_}, strike_, t_, ts_, {nMin_, nMax_}] /; Or[bUnd < bLow, bUnd > bHigh] :=
    doubleKnockOutOutsideBarrierOptionValue[cp, {sUnd, sYield, sVol}, {sUnd, sYield, sVol}, $unitCorrelation, intRate, {$zeroTime, t, $zeroBarrier, $infiniteBarrier}, strike, t, ts, {nMin, nMax}] - doubleKnockOutOutsideBarrierOptionValue[cp, {sUnd, sYield, sVol}, {bUnd, bYield, bVol}, corr, intRate, {t1, t2, bLow, bHigh}, strike, t, ts, {nMin, nMax}]

singleUpAndOutOutsideBarrierOptionValue[cp_, {sUnd_, sYield_, sVol_}, {bUnd_, bYield_, bVol_}, corr_, intRate_, {t1_, t2_, barrier_}, strike_, t_, ts_, {nMin_, nMax_}] /; barrier < bUnd :=
    0
singleUpAndOutOutsideBarrierOptionValue[cp_, {sUnd_, sYield_, sVol_}, {bUnd_, bYield_, bVol_}, corr_, intRate_, {t1_, t2_, barrier_}, strike_, t_, ts_, {nMin_, nMax_}] :=
    doubleKnockOutOutsideBarrierOptionValue[cp, {sUnd, sYield, sVol}, {bUnd, bYield, bVol}, corr, intRate, {t1, t2, $zeroBarrier, barrier}, strike, t, ts, {nMin, nMax}]

singleDownAndOutOutsideBarrierOptionValue[cp_, {sUnd_, sYield_, sVol_}, {bUnd_, bYield_, bVol_}, corr_, intRate_, {t1_, t2_, barrier_}, strike_, t_, ts_, {nMin_, nMax_}] /; bUnd < barrier :=
    0
singleDownAndOutOutsideBarrierOptionValue[cp_, {sUnd_, sYield_, sVol_}, {bUnd_, bYield_, bVol_}, corr_, intRate_, {t1_, t2_, barrier_}, strike_, t_, ts_, {nMin_, nMax_}] :=
    doubleKnockOutOutsideBarrierOptionValue[cp, {sUnd, sYield, sVol}, {bUnd, bYield, bVol}, corr, intRate, {t1, t2, barrier, $infiniteBarrier}, strike, t, ts, {nMin, nMax}]

singleDownAndInOutsideBarrierOptionValue[cp_, {sUnd_, sYield_, sVol_}, {bUnd_, bYield_, bVol_}, corr_, intRate_, {t1_, t2_, barrier_}, strike_, t_, ts_, {nMin_, nMax_}] /; bUnd < barrier :=
    0
singleDownAndInOutsideBarrierOptionValue[cp_, {sUnd_, sYield_, sVol_}, {bUnd_, bYield_, bVol_}, corr_, intRate_, {t1_, t2_, barrier_}, strike_, t_, ts_, {nMin_, nMax_}] :=
    doubleKnockOutOutsideBarrierOptionValue[cp, {sUnd, sYield, sVol}, {sUnd, sYield, sVol}, $unitCorrelation, intRate, {$zeroTime, t, $zeroBarrier, $infiniteBarrier}, strike, t, ts, {nMin, nMax}] - doubleKnockOutOutsideBarrierOptionValue[cp, {sUnd, sYield, sVol}, {bUnd, bYield, bVol}, corr, intRate, {t1, t2, barrier, $infiniteBarrier}, strike, t, ts, {nMin, nMax}]

singleUpAndInOutsideBarrierOptionValue[cp_, {sUnd_, sYield_, sVol_}, {bUnd_, bYield_, bVol_}, corr_, intRate_, {t1_, t2_, barrier_}, strike_, t_, ts_, {nMin_, nMax_}] /; barrier < bUnd :=
    0
singleUpAndInOutsideBarrierOptionValue[cp_, {sUnd_, sYield_, sVol_}, {bUnd_, bYield_, bVol_}, corr_, intRate_, {t1_, t2_, barrier_}, strike_, t_, ts_, {nMin_, nMax_}] :=
    doubleKnockOutOutsideBarrierOptionValue[cp, {sUnd, sYield, sVol}, {sUnd, sYield, sVol}, $unitCorrelation, intRate, {$zeroTime, t, $zeroBarrier, $infiniteBarrier}, strike, t, ts, {nMin, nMax}] - doubleKnockOutOutsideBarrierOptionValue[cp, {sUnd, sYield, sVol}, {bUnd, bYield, bVol}, corr, intRate, {t1, t2, $zeroBarrier, barrier}, strike, t, ts, {nMin, nMax}]

End[]

EndPackage[]

