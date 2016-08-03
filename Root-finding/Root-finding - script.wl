(* ::Package:: *)

Problem Set 2


Mateo Ochoa


Problem 3.11


Bisect[F_, xl_, xr_, eps_] := Module[{val, vam, xm,,xll,xrr},
xll=xl;
xrr=xr;
xm = .5(xll + xrr);
vam=F[xm];
While[Abs[vam]>eps,
xm = .5(xll + xrr);
val = F[xll];
vam = F[xm];
If[val*vam<0,
xrr=xm,
xll=xm];
];
Return[xm];]



BisectStart[F_,dx_,x_]:=Module[{N,interval},
N=1;
While[F[x+(N-1)dx]*F[x+N*dx]>0,
N=N+1];
interval={x+(N-1)dx,x+N*dx};
Return[interval];
]


Bisection[F_,deltax_,xstart_,eps_]:=Module[{interval,xm},
interval=BisectStart[F,deltax,xstart];
xm=Bisect[F,interval[[1]],interval[[2]],eps];
Return[xm]]


Problem 3.12


The Bisect function above modified to return the number of steps.


Bisect[F_, xl_, xr_, eps_] := Module[{N,val, vam, xm,,xll,xrr},
N=0;
xll=xl;
xrr=xr;
xm = .5(xll + xrr);
vam=F[xm];
While[Abs[vam]>eps,
xm = .5(xll + xrr);
val = F[xll];
vam = F[xm];
If[val*vam<0,
xrr=xm,
xll=xm];
N=N+1];
Return[N];
]


The Bisection function that will take the modified Bisect function as an input.


Bisection[F_,deltax_,xstart_,eps_]:=Module[{interval,xm},
interval=BisectStart[F,deltax,xstart];
xm=Bisect[F,interval[[1]],interval[[2]],eps];
Return[xm]]


F[x_]:=Sin[x^2-2*x-1]


Bisection[F,1,0,10^-10]


Newton[F_,DF_,x0_,eps_] := Module[{x, deltax,val},
x = x0;
val = F[x];
While[Abs[val] >= eps,
deltax = -F[x]/DF[x];
x = x + deltax;
val = F[x];
];
Return[x];
]


Newton[F,F',2,10^-10]


The Bisection function takes 29 steps to find the root, while the Newton value takes 6 steps.


Problem 3.14


G[x_]:=Tan[x]-(2*Sqrt[(100/x^2)-1])/(2-(100/x^2))


Plot[G[x],{x,0,12}]


Bisect[G,2.,4.,10^(-9)]


Bisect[G,5.,6.,10^(-9)]


Bisect[G,7.5,8.,10^(-9)]


Newton[G,G',2.,10^(-9)]


Newton[G,G',5.,10^(-9)]


Newton[G,G',7.5,10^(-9)]


The same results were obtained by using two different methods. 


Problem 3.18


RKODE4[G_, x0_,f0_, dx_, Ns_] := Module[{k1,k2,k3,k4,retvals,index,xn,fn},
retvals = Table[{0.0,0.0},{i,1,Ns}];
xn = x0;
fn = f0;
retvals[[1]] = {x0,f0};
For[index = 2, index <= Ns, index=index+1,
k1 = dx G[xn,fn];
k2 = dx G[xn + dx/2, fn + k1/2];
k3=dx G[xn+dx/2,fn+k2/2];
k4=dx G[xn+dx,fn+k3];
fn = fn + 1/3(k1/2 +k2+k3+k4/2);
xn = xn + dx;
retvals[[index]] = {xn,fn};
];
Return[retvals];
]


G[q_,f_]:={f[[2]],-(2/q)(f[[2]])+(-(2/q)-f[[3]])f[[1]],0}


F[Et_] := Module[{Nsteps,G,nvals,PsiE},
Nsteps = 20000;
G[x_,f_] := {f[[2]],-(2/x)(f[[2]])+(-(2/x)-f[[3]])f[[1]],0};
nvals = RKODE4[G, 20.0/(Nsteps-1), {0.0,1.0,Et},25.0/(Nsteps-1), Nsteps];
PsiE = nvals[[Length[nvals],2,1]];
Return[PsiE];
]


25/19999.


The values for the first three bound states for hydrogen and their accepted values 1/n^2 will be:


Bisect[F,-1.5,-0.9,10^(-9)]


The accepted value is -1.


Bisect[F,-0.9913449126383057`,-0,10^(-9)]


The accepted value is -0.25.


Bisect[F,-0.23,-0,10^(-9)]


The accepted value is -0.11.


Changing the value for Subscript[q, infinity] to 25 and Nsteps to 20000, the first three bound states will be:


Bisect[F,-1.5,-0.9,10^(-9)]


Bisect[F,-0.8,-0,10^(-9)]


Bisect[F,-0.23,-0,10^(-9)]


We can see that the latter approximations are more accurate than the former ones. 
