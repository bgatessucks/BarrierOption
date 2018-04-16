(* Wolfram Language Package *)

(* Created by the Wolfram Workbench 10-Apr-2015 *)

BeginPackage["BarrierOption`"]
(* Exported symbols added here with SymbolName::usage *) 


bp::usage="Helper function"
bpp::usage="Helper function"

(*defgGeneric::usahe=""*)
(*phi::usage="Utility function"*)

barrierOption::usage="Returns the value of a barrier option"
doubleKnockOutOutsideBarrierOptionValue::usage="Returns the value of a double knock-out outside barrier option."

Begin["`Private`"]
(* Implementation of the package *)

daysYear = 365.242;

bp[low_, high_, n_] = 2 (high - low) n;
bpp[low_, high_, n_] = 2 high -  bp[low, high, n];

defgGeneric[nontterm_, tterm_, t_, lastterm_, den_, additionalterm_] = (nontterm + tterm t + lastterm) / den + additionalterm;

phi = Function[{mean, cov}, CDF[MultinormalDistribution[mean, cov], #] &];
SetAttributes[phi, NumericFunction];

(* CDF[MultinormalDostribution[, ], ] *)

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

doubleKnockOutOutsideBarrierOptionValue[cp_, {sUnd_, sYield_, sVol_}, {bUnd_, bYield_, bVol_}, corr_, intRate_, {t1_, t2_, bLow_, bHigh_}, strike_, t_, ts_, {nMin_, nMax_}] /; Or[bUnd < bLow, bUnd > bHigh] :=
    0


doubleKnockOutOutsideBarrierOptionValue[cp_, {sUnd_, sYield_, sVol_}, {bUnd_, bYield_, bVol_}, corr_, intRate_, {t1_, t2_, bLow_, bHigh_}, strike_, t_, ts_, {nMin_, nMax_}] :=
    Module[ { t1Local, n, 
            a, bH, bL, mu1, mu2, 
            gMat, d1p, d2p, d1pp, d2pp, e1p, e2p, e1pp, e2pp, f1p, f2p, f1pp, f2pp, g1, g2, h1, h2,
            wsn, wkn, ws, wk, 
            value
        },
        t1Local = Max[t1, $MachineEpsilon]; (* Takes care of t1=0 in g1, g2, h1, h2 *)
        
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
        ws = Sum[Evaluate[wsn], {n, nMin, nMax}];
        wk = Sum[Evaluate[wkn], {n, nMin, nMax}];
        value = Exp[-intRate (ts - t)] cp (sUnd Exp[-sYield t] ws - strike Exp[-intRate t] wk)
    ]

End[]

EndPackage[]

