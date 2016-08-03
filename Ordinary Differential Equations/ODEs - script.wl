(* ::Package:: *)

(* ::Title:: *)
(*Mateo Ochoa*)


(* ::Section:: *)
(*Lab 1*)


(* ::Subsection:: *)
(*Problem 2.15*)


(* ::Subsubsection:: *)
(*Part a*)


(* ::Input:: *)
(*Verlet[a_, x0_, v0_, dt_, Ns_] := Module[{xout,xnext, xcur,xprev,t,index},*)
(*xout = Table[0.0,{j,1,Ns}];*)
(*xcur = x0;*)
(*xprev = x0 - dt v0 +(1/2) a[x0] dt^2;*)
(*t = 0.0;*)
(*xout[[1]] = {t,x0};*)
(*For[index = 2, index <= Ns, index=index+1,*)
(*xnext = 2.0 xcur - xprev + dt^2 a[xcur];*)
(*t = t + dt;*)
(*xout[[index]] = {t, xnext};*)
(*xprev = xcur;*)
(*xcur = xnext;*)
(*];*)
(*Return[xout];*)
(*]*)


(* ::Input:: *)
(*a[x_]:=-(2Pi)^2*x*)


(* ::Input:: *)
(*listx=Verlet[a,1,0,20/(999),1000];*)


(* ::Input:: *)
(*v=Table[{listx[[i+1,1]],(listx[[i+2,2]]-listx[[i,2]])/(2*20/(999))},*)
(*{i,1,Length[listx]-2}];*)


(* ::Input:: *)
(*energy=Table[{listx[[i+1,1]],((v[[i,2]]^2)/2)+((2Pi)^2*(listx[[i+1,2]])^2)/2},*)
(*{i,1,Length[v]-1}];*)


(* ::Input:: *)
(*ListPlot[energy]*)


(* ::Input:: *)
(*Sum[energy[[i,2]],{i,Length[energy]}]/Length[energy]*)


(* ::Input:: *)
(*Max[Table[energy[[i,2]],{i,1,Length[energy]}]]-Min[Table[energy[[i,2]],{i,1,Length[energy]}]]*)


(* ::Subsubsection:: *)
(*Part b*)


(* ::Input:: *)
(*RKODE[G_, t0_,x0_, dt_, Ns_] := Module[{k1,k2,retvals,index,tn,xn},*)
(*retvals = Table[{0.0,0.0},{i,1,Ns}];*)
(*tn = t0;*)
(*xn = x0;*)
(*retvals[[1]] = {t0,x0};*)
(*For[index = 2, index <= Ns, index=index+1,*)
(*k1 = dt G[tn,xn];*)
(*k2 = dt G[tn + dt/2, xn + k1/2];*)
(*xn = xn + k2;*)
(*tn = tn + dt;*)
(*retvals[[index]] = {tn,xn};*)
(*];*)
(*Return[retvals];*)
(*]*)


(* ::Input:: *)
(*G[t_,f_] :={f[[2]],-(2Pi)^2*f[[1]]}*)


(* ::Input:: *)
(*v2=RKODE[G,0.,{1.,0.},20/(999),1000];*)


(* ::Input:: *)
(*xtab = Table[{v2[[j,1]],v2[[j,2,1]]},{j,1,Length[v2]}];*)


(* ::Input:: *)
(*ListPlot[xtab]*)


(* ::Input:: *)
(*vtab = Table[{v2[[j,1]],v2[[j,2,2]]},{j,1,Length[v2]}];*)


(* ::Input:: *)
(*energy2=Table[{xtab[[i,1]],((vtab[[i,2]]^2)/2)+((2Pi)^2*(xtab[[i,2]])^2)/2},*)
(*{i,1,Length[v]-1}];*)


(* ::Input:: *)
(*ListPlot[energy2]*)


(* ::Subsection:: *)
(*Problem 2.17*)


(* ::Input:: *)
(*RKODE[G_, x0_,f0_, dx_, Ns_] := Module[{k1,k2,retvals,index,xn,fn},*)
(*retvals = Table[{0.0,0.0},{i,1,Ns}];*)
(*xn = x0;*)
(*fn = f0;*)
(*retvals[[1]] = {x0,f0};*)
(*For[index = 2, index <= Ns, index=index+1,*)
(*k1 = dx G[xn,fn];*)
(*k2 = dx G[xn + dx/2,fn+k1/2];*)
(*fn = fn + k2;*)
(*xn = xn + dx;*)
(*retvals[[index]] = {xn,fn};*)
(*];*)
(*Return[retvals];*)
(*]*)


(* ::Input:: *)
(*RKODE4[G_, x0_,f0_, dx_, Ns_] := Module[{k1,k2,k3,k4,retvals,index,xn,fn},*)
(*retvals = Table[{0.0,0.0},{i,1,Ns}];*)
(*xn = x0;*)
(*fn = f0;*)
(*retvals[[1]] = {x0,f0};*)
(*For[index = 2, index <= Ns, index=index+1,*)
(*k1 = dx G[xn,fn];*)
(*k2 = dx G[xn + dx/2, fn + k1/2];*)
(*k3=dx G[xn+dx/2,fn+k2/2];*)
(*k4=dx G[xn+dx,fn+k3];*)
(*fn = fn + 1/3(k1/2 +k2+k3+k4/2);*)
(*xn = xn + dx;*)
(*retvals[[index]] = {xn,fn};*)
(*];*)
(*Return[retvals];*)
(*]*)


(* ::Subsubsection:: *)
(*Part a*)


(* ::Input:: *)
(*g1[x_,f_]:=x^10-5x^2*)


(* ::Input:: *)
(*eq1=RKODE[g1,0,1,0.001,1/0.001];*)


(* ::Input:: *)
(*s1[x_] := x^11/11 - 5/3 x^3+1*)


(* ::Input:: *)
(*res1 = Table[{eq1[[j,1]],s1[eq1[[j,1]]] - eq1[[j,2]]},{j,1,Length[eq1]}];*)


(* ::Input:: *)
(*ListPlot[res1]*)


(* ::Input:: *)
(*eq11=RKODE4[g1,0,1,0.001,1/0.001];*)


(* ::Input:: *)
(*res11 = Table[{eq11[[j,1]],s1[eq11[[j,1]]] - eq11[[j,2]]},{j,1,Length[eq11]}];*)


(* ::Input:: *)
(*ListPlot[res11]*)


(* ::Subsubsection:: *)
(*Part b*)


(* ::Input:: *)
(*g2[x_,f_]:=-x f*)


(* ::Input:: *)
(*eq2=RKODE[g2,0,1,0.001,1/0.001];*)


(* ::Input:: *)
(*s2[x_]:=Exp[(-x^2/2)]*)


(* ::Input:: *)
(*res2 = Table[{eq2[[j,1]],s2[eq2[[j,1]]] - eq2[[j,2]]},{j,1,Length[eq2]}];*)


(* ::Input:: *)
(*ListPlot[res2]*)


(* ::Input:: *)
(*eq22=RKODE4[g2,0,1,0.001,1/0.001];*)


(* ::Input:: *)
(*res22 = Table[{eq22[[j,1]],s2[eq22[[j,1]]] - eq22[[j,2]]},{j,1,Length[eq22]}];*)


(* ::Input:: *)
(*ListPlot[res22]*)


(* ::Subsubsection:: *)
(*Part c*)


(* ::Input:: *)
(*g3[x_,f_]:={f[[2]],-29 f[[1]] -4 f[[2]]}*)


(* ::Input:: *)
(*eq3=RKODE[g3,0,{1,0},0.001,1/0.001];*)


(* ::Input:: *)
(*s3[x_]:=(1/2 -1/5 I) Exp[(-2+5 I)x]+(1/2 +1/5 I) Exp[(-2-5 I)x]*)


(* ::Input:: *)
(*res3 = Table[{eq3[[j,1]],Chop[s3[eq3[[j,1]]] ]- eq3[[j,2,1]]},{j,1,Length[eq3]}];*)


(* ::Input:: *)
(*ListPlot[res3]*)


(* ::Input:: *)
(*eq33=RKODE4[g3,0,{1,0},0.001,1/0.001];*)


(* ::Input:: *)
(*res33= Table[{eq33[[j,1]],Chop[s3[eq33[[j,1]]] ]- eq33[[j,2,1]]},{j,1,Length[eq33]}];*)


(* ::Input:: *)
(*ListPlot[res33]*)


(* ::Input:: *)
(*eq3[[1+.95/0.001,2,1]]*)


(* ::Input:: *)
(*eq33[[1+.95/0.001,2,1]]*)


(* ::Subsection:: *)
(*Problem 2.19*)


(* ::Subsubsection:: *)
(*Part a*)


(* ::Input:: *)
(*c=3*10^8;*)


(* ::Input:: *)
(*relspring[x_,f_]:={f[[2]],-12 f[[1]] (1-f[[2]]^2/c^2)^(3/2)}*)


(* ::Input:: *)
(*rs=RKODE4[relspring,0,{1,0},0.05,10/0.05];*)


(* ::Input:: *)
(*xrs=Table[{rs[[i,1]],rs[[i,2,1]]},{i,1,Length[rs]}];*)


(* ::Input:: *)
(*xns[x_]:=Cos[Sqrt[12] x]*)


(* ::Input:: *)
(*Show[Plot[xns[x],{x,0,10}],ListPlot[xrs]]*)


(* ::Input:: *)
(*vrs=Table[{rs[[i,1]],rs[[i,2,2]]},{i,1,Length[rs]}];*)


(* ::Input:: *)
(*vns[x_]:=-Sqrt[12] Sin[Sqrt[12] x]*)


(* ::Input:: *)
(*Show[Plot[vns[x],{x,0,10}],ListPlot[vrs]]*)


(* ::Text:: *)
(*Yes, the relativistic solution approximates the non-relativistic one.*)


(* ::Subsubsection:: *)
(*Part b*)


(* ::Input:: *)
(*rs2=RKODE4[relspring,0,{c/Sqrt[12],0},0.05,10/0.05];*)


(* ::Input:: *)
(*xrs2=Table[{rs2[[i,1]],rs2[[i,2,1]]},{i,1,Length[rs2]}];*)


(* ::Input:: *)
(*xns2[x_]:=(c/Sqrt[12])Cos[Sqrt[12] x]*)


(* ::Input:: *)
(*Show[Plot[xns2[x],{x,0,10}],ListPlot[xrs2]]*)


(* ::Text:: *)
(*The solid line is the newtonian approach, whereas the dotted line is the relativistic approach.*)


(* ::Input:: *)
(*vrs2=Table[{rs2[[i,1]],rs2[[i,2,2]]},{i,1,Length[rs2]}];*)


(* ::Input:: *)
(*vns2[x_]:=-(c/Sqrt[12])Sqrt[12] Sin[Sqrt[12] x]*)


(* ::Input:: *)
(*Show[Plot[vns2[x],{x,0,10}],ListPlot[vrs2]]*)


(* ::Text:: *)
(*The solid line is the newtonian approach, whereas the dotted line is the relativistic approach.*)


(* ::Input:: *)
(*Max[Table[{vrs2[[i,2]]},{i,1,Length[xrs2]}]]*)


(* ::Subsubsection:: *)
(*Part c*)


(* ::Input:: *)
(*rs3=RKODE4[relspring,0,{5*c/Sqrt[12],0},0.05,10/0.05];*)


(* ::Input:: *)
(*xrs3=Table[{rs3[[i,1]],rs3[[i,2,1]]},{i,1,Length[rs3]}];*)


(* ::Input:: *)
(*xns3[x_]:=5*(c/Sqrt[12])Cos[Sqrt[12] x]*)


(* ::Input:: *)
(*Show[Plot[xns3[x],{x,0,10}],ListPlot[xrs3]]*)


(* ::Text:: *)
(*The solid line is the newtonian approach, whereas the dotted line is the relativistic approach.*)


(* ::Input:: *)
(*vrs3=Table[{rs3[[i,1]],rs3[[i,2,2]]},{i,1,Length[rs3]}];*)


(* ::Input:: *)
(*vns3[x_]:=-(c/Sqrt[12])Sqrt[12] Sin[Sqrt[12] x]*)


(* ::Input:: *)
(*Show[Plot[vns3[x],{x,0,10}],ListPlot[vrs3]]*)


(* ::Text:: *)
(*The solid line is the newtonian approach, whereas the dotted line is the relativistic approach.*)


(* ::Input:: *)
(*Max[Table[{vrs3[[i,2]]},{i,1,Length[xrs3]}]]*)


(* ::Text:: *)
(*From the position graph we see what resembles a triangle wave. From the velocity graph we see that the velocity does not change much in magnitude, but changes its direction periodically. The mass travels from one extreme to the other at almost constant speed; it only slows down when its displacement is very close to the amplitude. *)


(* ::Subsection:: *)
(*Problem 2.20*)


(* ::Subsubsection:: *)
(*Part a*)


(* ::Input:: *)
(*motion1[t_,f_]:={f[[2]],0,f[[4]],-9.8}*)


(* ::Input:: *)
(*motion11=RKODE4[motion1,0,{0,1000*Cos[Pi/6],0,1000*Sin[Pi/6]},0.04,2552];*)
(*coord=Table[{motion11[[i,2,1]],motion11[[i,2,3]]},{i,0,Length[motion11]}];*)
(*ListPlot[coord]*)


(* ::Text:: *)
(*To obtain the appropiate step size, the total time of flight was obtained:*)


(* ::Input:: *)
(*(2000*Sin[Pi/6])/9.8*)


(* ::Text:: *)
(*The number of steps to represent the whole trajectory was this number +1:*)


(* ::Input:: *)
(*102.0408163265306`/0.04*)


(* ::Text:: *)
(*The last element has its time component with the value of*)


(* ::Input:: *)
(*motion11[[2552,1]]*)


(* ::Text:: *)
(*which is accurate to the first two decimals.*)


(* ::Subsubsection:: *)
(*Part b*)


(* ::Input:: *)
(*motion2[t_,f_]:={f[[2]],-(5*10^(-5))*Sqrt[f[[2]]^2+f[[4]]^2]*f[[2]],f[[4]],-9.8-(5*10^(-5))*Sqrt[f[[2]]^2+f[[4]]^2]*f[[4]]}*)


(* ::Input:: *)
(*motion21=RKODE4[motion2,0,{0,1000*Cos[Pi/6],0,1000*Sin[Pi/6]},0.04,1800];*)
(*coord=Table[{motion21[[i,2,1]],motion21[[i,2,3]]},{i,0,Length[motion21]}];*)
(*ListPlot[coord]*)


(* ::Subsubsection:: *)
(*Part c*)


(* ::Input:: *)
(*motion3[t_,f_]:={f[[2]],-(5*10^(-5))*Sqrt[f[[2]]^2+f[[4]]^2]*f[[2]]-.5*f[[4]],f[[4]],-9.8-(5*10^(-5))*Sqrt[f[[2]]^2+f[[4]]^2]*f[[4]]+.5*f[[2]]}*)


(* ::Input:: *)
(*motion31=RKODE4[motion3,0,{0,45*Cos[Pi/18],0,45*Sin[Pi/18]},0.04,290];*)
(*coord=Table[{motion31[[i,2,1]],motion31[[i,2,3]]},{i,0,Length[motion31]}];*)
(*ListPlot[coord]*)
