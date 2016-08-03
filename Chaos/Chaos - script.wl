(* ::Package:: *)

Lab 3


Mateo Ochoa


Problem 13.11


recursion[G_,uo_,N_,a_]:=Module[{n,u},
n=1;
u=uo;
While[n<=N,
u=G[u];
n=n+1;
];
Return[u];
]


g[x_,a_]:=3x-x^3


Recursion[G_,uo_,a_,N_,de_,ei_,ef_]:=Module[{n,nn,u,e,list,plot},
nn=1;
list=Table[{0,0},{i,(ef-ei)/de +1}];
e=ei;
While[e<=ef,
n=1;
u=uo+e;
While[n<=N,
u=G[u,a];
n=n+1;
];
list[[nn]]={e,u};
e=e+de;
nn=nn+1;
];
plot=ListPlot[list];
Return[plot];
]


Recursion[g,Sqrt[2.],100,10^(-11),10^(-10),10^(-8)]


We can see that Sqrt[2] is not a stable fixed point.


Problem 13.15


The values of u of a function G will be given by:


valu[G_,uo_,a_,N_]:=Module[{list,n,u},
n=2;
list=Table[0,{i,1,N}];
list[[1]]=uo;
u=uo;
While[n<=N,
u=G[u,a];
list[[n]]=u;
n=n+1;
];
Return[list]
]


The fixed points will be taken from the last m u values.


fixedp[values_,m_,e_]:=Module[{last,fixed},
last=Take[values,-m];
fixed=Union[last, SameTest-> (Abs[#1-#2] <= e&)];
Return[fixed];
]


For a range of values of a, the fixed points will be given:


avsu[G_,uo_,ai_,af_,da_,N_,m_,e_]:=Module[{values,fxdp,a,list,n,plot},
list=Table[0,{i,1,(af-ai)/da +1}];
a=ai;
n=1;
While[a<=af,
values=valu[G,uo,a,N];
fxdp=fixedp[values,m,e];
list[[n]]={a,fxdp};
a=a+da;
n=n+1;

];
Return[list];
]


PlotFurcate[avals_] := Module[{ret,pts,index,entry,jndex},
pts = {PointSize[0.001]};
For[index = 1, index <= Length[avals], index=index+1,
entry = avals[[index]];
pts = Join[pts,Table[Point[{entry[[1]], entry[[2,j]]}],{j,1,Length[entry[[2]]]}]];
];
ret = Show[Graphics[{PointSize[.01/2],pts}],Axes->True,AxesLabel->{"a","u*"},AxesOrigin->{1,0}];
Return[ret];
]


Defining the function f(x):


f[x_,a_]:=a Sin[Pi x]


avals=avsu[f,Pi/4.,0.7,0.86,0.0001,200,20,10^(-4)];


The bifurcation diagram of the function f will be:


PlotFurcate[avals]


The two closest a values to Subscript[a, 1] that only have two fixed points are:


{0.7112999999999987`,{0.6417918702754514`,0.6418921089823065`}}
{0.7223999999999975`,{0.6155933344529679`,0.6752173884264112`}}


Taking the average, the value of Subscript[a, 1] will be:


a1=(0.7112999999999987+0.7223999999999975)/2


Similarly, for Subscript[a, 2] and Subscript[a, 3]:


a2=(0.8302999999999856`+ 0.8303999999999856`)/2


a3=0.8577999999999826


The approximation to the Feigenbaum number will be:


\[Delta]=(a2-a1)/(a3-a2)


Problem 13.16


g[x_,a_]:=a x- a x^2 


This function gives the coefficient for a given a:


Lyap[G_,uo_,a_,N_]:=Module[{n,lambda,,nn,u,e,list,plot,\[Lambda]},
nn=2;
u=uo;
list=Table[{0,0},{i,1,N+2}];
lambda=Log[Abs[a-2 a u]];
list[[1]]={uo,lambda};
n=0;
While[n<=N,
u=G[u,a];
n=n+1;
lambda=lambda+Log[Abs[a - 2*a*u]];
list[[nn]]={u,lambda};
nn=nn+1
];
\[Lambda]=list[[N+2,2]]/(N+2);
Return[\[Lambda]];
]


Testing the function for a know value (a = 2.5,  given in the book):


Lyap[g,0.1,2.5,1000]


This function gives the values of the Lyapunov exponent for a range of values of a:


Lyapvsa[G_,uo_,ai_,af_,da_,N_]:=Module[{list,a,coeff,n,plot},
list=Table[{0,0},{i,1,(af-ai)/da +1}];
a=ai;
n=1;
While[a<=af,
coeff=Lyap[G,uo,a,N];
list[[n]]={a,coeff};
a=a+da;
n=n+1;
];
plot=ListPlot[list,AxesLabel->{"a",\[Lambda]}];
Return[plot];
]


Lyapvsa[g,0.1,2.9,4,0.01,1000]


When the exponent becomes positive, there are no stable points for those values of a. At the bifurcation points, the exponent decreases significantly. This happens because at the bifurcation points the difference between the values of the fixed points decreases. In other words, the "slope" of the bifurcation graph approaches zero. When the exponent has negative spikes, it means that fixed points exist for those values of a.


Problem 13.19


Multibrot[q_,xi_,xf_,yi_,yf_,N_,Npts_]:=Module[{plot,list,zo,z,n,X,Y,deltax,deltay,c},
deltax=(xf-xi)/Npts;
deltay=(yf-yi)/Npts;
list={};
For[X = xi, X <= xf, X = X + deltax,
For[Y = yi, Y <= yf, Y = Y + deltay,
zo=X + I Y;
z = zo;
n = 1;
While[n<=N,
z=z^(q) +zo;
n=n+1;
];
If[Abs[z] < 10,list = Join[list,{{X,Y}}];
];
];
];
Return[list];
]


test = Multibrot[3,-1.,1.,-3.,3.,10,1000];


Show[Graphics[Table[Point[test[[j]]],{j,1,Length[test]}]]]
