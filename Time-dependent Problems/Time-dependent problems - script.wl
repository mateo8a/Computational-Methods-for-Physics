(* ::Package:: *)

Lab 5


Mateo Ochoa


Problem 5.11


Part a


This function follows the Richtmyer two-step method in full non-linear form.


uhalf[u0_,v_,dt_,Nt_,dx_,bndtype_]:=Module[{Nx,uup,u,upos,i,j,uplus12,uminus12,uplus1},
i=1;
j=2;
u=u0;
Nx=Length[u0];
upos=Table[0.0,{i,1,Nx}];
uup=Table[0.0,{j,1,Nt}];

While[i<=Nt,
While[j<Nx,
uplus12=(u[[j]]+u[[j+1]])/2 -(dt/(2 dx)) (u[[j+1]] v[dx(j+1),dt i]-u[[j]]v[dx j,dt i]);
uminus12=(u[[j-1]]+u[[j]])/2 -(dt/(2 dx)) (u[[j]] v[dx j,dt i]-u[[j-1]]v[dx(j-1),dt i]);
uplus1=u[[j]]-(dt/dx)(uplus12 v[dx(j+1/2),dt(i+1/2)] - uminus12 v[dx(j-1/2),dt(i+1/2)]);
upos[[j]]=uplus1;
j=j+1];

If[bndtype=="periodic",
upos[[1]] = upos[[Nx-1]];
upos[[Nx]] = upos[[2]];
];
If[bndtype =="diffuse",
upos[[1]] = upos[[2]];
upos[[Nx]] = upos[[Nx-1]];
];
If[bndtype =="zero",
upos[[1]] = 0.0;
upos[[Nx]] = 0.0;
];

uup[[i]]={upos,i dt};
i=i+1;
j=2;
u=upos];

Return[uup];
]


Let the conditions for the problem be:


dx=0.01;
v[x_,t_]:=50;
dt=dx/v[0.,0.];
Nt=200;
f[x_]:=Sin[Pi x];
u0=Table[f[i],{i,0,1+dx,dx}];


Then, a table with the plots of f[x] for every t = j dt will be created. With this in place, we animate the list of plots. 


wave=uhalf[u0,v,dt,Nt,dx,"periodic"];
mov = Table[ListPlot[wave[[j,1]]],{j,1,Length[wave]}];
ListAnimate[mov]


Part b


The function from part a will be used. This time, however, the spacing in the spatial grid, dx, will be different; dx will be 1.1 times greater than in part a.


dx=0.01;
v[x_,t_]:=50;
dt=1.1 dx/v[0.,0.];
Nt=200;
f[x_]:=Sin[Pi x];
u0=Table[f[i],{i,0,1+dx,dx}];


wave=uhalf[u0,v,dt,Nt,dx,"periodic"];
mov = Table[ListPlot[wave[[j,1]]],{j,1,Length[wave]}];
ListAnimate[mov]


Problem 5.14


Now we'll solve the traffic flow problem.


dx=0.01;
v[x_,t_]:=65;
dt=dx/v[0.,0.];
Nt=200;
f[x_]:=If[x<499dx,0,132];
u0=Table[f[i],{i,0,10,dx}];


traffic=uhalf[u0,v,dt,Nt,dx,"diffuse"];


mov = Table[ListPlot[traffic[[j,1]]],{j,1,Length[traffic]}];
ListAnimate[mov]


From the animation we can approximate the speed. 


(700-500)/(Nt dt)


The theoretical value for the speed of the shock is


F[t_]:=(1-(132/264)) v[0.,0.] t


F[Nt]


as expected.


We will solve the problem using Lax-Friedrichs an compare with the result obtained using Lax-Wendroff.


LFPDE[u0_, v_, dt_,Nt_,dx_, bndtype_] := Module[{retu, u,uup, Nx,i,x,j},
u = u0;
Nx = Length[u0];
uup = Table[0,{k,1,Nx}];
retu = Table[0,{k,1,Nt}];
For[i = 1, i <= Nt, i = i+1,
x= dx;
For[j = 2, j < Nx, j = j+1,
uup[[j]] = (1/2) (u[[j+1]] + u[[j-1]]) 
- dt/(2 dx) ( u[[j+1]] v[x+dx,u[[j+1]]] - u[[j-1]] v[x-dx,u[[j-1]]]);
x = x + dx;
];
If[bndtype=="periodic",
uup[[1]] = uup[[Nx-1]];
uup[[Nx]] = uup[[2]];
];
If[bndtype =="diffuse",
uup[[1]] = uup[[2]];
uup[[Nx]] = uup[[Nx-1]];
];
If[bndtype =="zero",
uup[[1]] = 0.0;
uup[[Nx]] = 0.0;
];
u = uup;
retu[[i]] = uup;
];
Return[retu];
]



traffic=LFPDE[u0,v,dt,Nt,dx,"diffuse"];
mov = Table[ListPlot[traffic[[j]]],{j,1,Length[traffic]}];
ListAnimate[mov]


The speed using Lax-Friedrichs is:


(700-500)/(Nt dt)


The value obtained is the same value as the one obtained using Lax-Wendroff.


Problem 5.15


The explicit Euler method is implemented below.


This will define the matrix that will operate on Subscript[u, j]^nfor j \[Epsilon] [1,Nq].


ExplicitE[u0_, v_, dt_,dq_] := Module[{g,Nq,retmat,indeq,indet,i,j,k},
Nq=Length[u0];

retmat = Table[0.0,{i,1,Nq },{j,1,Nq }];

retmat[[1,1]]=-2.0 dt I/(dq^2)-I dt v[dq indeq-1/2,indet dt]+1;
retmat[[1,2]] += I dt/dq^2;

For[indeq=2,indeq <= Nq-1, indeq = indeq+1,
i=indeq-1;
j=indeq;
k=indeq+1;
retmat[[j,i]] += I dt/dq^2;
retmat[[j,j]] += -2.0 dt I/(dq^2)-I dt v[dq indeq-1/2,indet dt]+1;
retmat[[j,k]] += I dt/dq^2;
];

retmat[[Nq,Nq-1]] += I dt/dq^2;
retmat[[Nq,Nq]] += -2.0 dt I/(dq^2)-I dt v[dq indeq-1/2,indet dt]+1;
Return[retmat];
]


The initial values of the problem are:


dq=1/201;
v[q_,s_]:=If[q<=1/8,0,5000];
ds=dq/2000;
Ns=800;
EE=4500;
p0=100;
f[q_]:= (2/Pi)^(1/4) (EE-(p0^2)/4)^(1/4) Exp[-(EE-(p0^2)/4)q^2] Exp[I (p0/2) q];
u0=Table[f[i],{i,-1/2+dq,1/2-dq,dq}];


This function will give a list of u^n for n \[Epsilon] [1,Ns]. Note that this will be a list where the elements (u for every n) are lists themselves (all the positions points j of u at n). 


ExplicitEtime[Ns_]:=Module[{Explicitu,TimeEv,n,u,Operator},
n=1;
TimeEv=Table[0.0,{i,1,Ns}];
u=u0;
Operator=ExplicitE[u,v,ds,dq];
While[n<= Ns,
Explicitu=Operator.u;
TimeEv[[n]]=Explicitu;
u=Explicitu;
n=n+1;
];
Return[TimeEv]
]


The values of u are complex numbers. In order to visualize the evolution of the wave, we compute the value of u*u to obtain real values. 


psi=ExplicitEtime[Ns];
psiabs=(Abs[psi])^2;
mov = Table[ListLinePlot[psiabs[[j]],PlotRange->{0,40}],{j,1,Length[psiabs]}];
ListAnimate[mov]


Problem 5,16


Now, the implicit Euler method will be implemented following the same steps as in Problem 5.15.


The matrix will be given by:


ImplicitE[u0_, v_, dt_,dq_] := Module[{g,Nq,retmat,indeq,indet,i,j,k},
Nq=Length[u0];

retmat = Table[0.0,{i,1,Nq },{j,1,Nq }];

retmat[[1,1]]=2.0 dt I/(dq^2)+I dt v[dq indeq-1/2,indet dt]+1;
retmat[[1,2]] += -I dt/dq^2;

For[indeq=2,indeq <= Nq-1, indeq = indeq+1,
i=indeq-1;
j=indeq;
k=indeq+1;
retmat[[j,i]] += -I dt/dq^2;
retmat[[j,j]] += 2.0 dt I/(dq^2)+I dt v[dq indeq-1/2,indet dt]+1;
retmat[[j,k]] += -I dt/dq^2;
];

retmat[[Nq,Nq-1]] += -I dt/dq^2;
retmat[[Nq,Nq]] += 2.0 dt I/(dq^2)+I dt v[dq indeq-1/2,indet dt]+1;
Return[retmat];
]


ImplicitEtime[Ns_]:=Module[{Implicitu,TimeEv,n,u,Operator},
n=1;
TimeEv=Table[0.0,{i,1,Ns}];
u=u0;
Operator=ImplicitE[u,v,ds,dq];
While[n<= Ns,
Implicitu=LinearSolve[Operator,u];
TimeEv[[n]]=Implicitu;
u=Implicitu;
n=n+1;
];
Return[TimeEv]
]


psi=ImplicitEtime[Ns];
psiabs=(Abs[psi])^2;
mov = Table[ListLinePlot[psiabs[[j]],PlotRange->{0,40}],{j,1,Length[psiabs]}];
ListAnimate[mov]


Problem 5.17


CN1[u0_, v_, dt_,dq_] := Module[{g,Nq,retmat,indeq,indet,i,j,k},
Nq=Length[u0];

retmat = Table[0.0,{i,1,Nq },{j,1,Nq }];

retmat[[1,1]]=-2.0 dt I/(2 dq^2)-I dt v[dq indeq-1/2,indet dt]/2+1;
retmat[[1,2]] += I dt/(2 dq^2);

For[indeq=2,indeq <= Nq-1, indeq = indeq+1,
i=indeq-1;
j=indeq;
k=indeq+1;
retmat[[j,i]] += I dt/(2 dq^2);
retmat[[j,j]] += -2.0 dt I/(2 dq^2)-I dt v[dq indeq-1/2,indet dt]/2+1;
retmat[[j,k]] += I dt/(2 dq^2);
];

retmat[[Nq,Nq-1]] += I dt/(2 dq^2);
retmat[[Nq,Nq]] += -2.0 dt I/(2 dq^2)-I dt v[dq indeq-1/2,indet dt]/2+1;
Return[retmat];
]


CN2[u0_, v_, dt_,dq_] := Module[{g,Nq,retmat,indeq,indet,i,j,k},
Nq=Length[u0];

retmat = Table[0.0,{i,1,Nq },{j,1,Nq }];

retmat[[1,1]]=+2.0 dt I/(2 dq^2)+I dt v[dq indeq-1/2,indet dt]/2+1;
retmat[[1,2]] +=- I dt/(2 dq^2);

For[indeq=2,indeq <= Nq-1, indeq = indeq+1,
i=indeq-1;
j=indeq;
k=indeq+1;
retmat[[j,i]] +=- I dt/(2 dq^2);
retmat[[j,j]] += 2.0 dt I/(2 dq^2)+I dt v[dq indeq-1/2,indet dt]/2+1;
retmat[[j,k]] += -I dt/(2 dq^2);
];

retmat[[Nq,Nq-1]] +=- I dt/(2 dq^2);
retmat[[Nq,Nq]] += 2.0 dt I/(2 dq^2)+I dt v[dq indeq-1/2,indet dt]/2+1;
Return[retmat];
]


CNtime[Ns_]:=Module[{Implicitu,TimeEv,n,u,Operator1,Operator2},
n=1;
TimeEv=Table[0.0,{i,1,Ns}];
u=u0;
Operator1=CN1[u,v,ds,dq];
Operator2=CN2[u,v,ds,dq];
While[n<= Ns,
Implicitu=LinearSolve[Operator2,Operator1.u];
TimeEv[[n]]=Implicitu;
u=Implicitu;
n=n+1;
];
Return[TimeEv]
]


psi=CNtime[Ns];
psiabs=(Abs[psi])^2;
mov = Table[ListLinePlot[psiabs[[j]],PlotRange->{0,40}],{j,1,Length[psiabs]}];
ListAnimate[mov]


Problem 5.18


Ns=500;


The norm values for the Explicit Euler method will be:


psiabs=(Abs[ExplicitEtime[Ns]])^2;
ListPlot[Table[Sum[psiabs[[j,i]]*dq,{i,200}],{j,1,Ns}],AxesLabel->{"Ns","norm value"}]


Then, for the Implicit Euler method:


psiabs=(Abs[ImplicitEtime[Ns]])^2;
ListPlot[Table[Sum[psiabs[[j,i]]*dq,{i,200}],{j,1,Ns}],AxesLabel->{"Ns","norm value"}]


Finally, for the Crank-Nicolson method:


psiabs=(Abs[CNtime[Ns]])^2;
ListPlot[Table[Sum[psiabs[[j,i]]*dq,{i,200}],{j,1,Ns}],AxesLabel->{"Ns","norm value"}]
