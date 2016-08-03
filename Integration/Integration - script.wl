(* ::Package:: *)

Lab 6


Mateo Ochoa


Problem 6.10


The following function is based off the "Lagrange's interpolating formula", with n=3.


LagrangeIF[f_,xi_,xf_,Ns_]:=Module[{w0,w1,w2,w3,dx,J,i},
dx=(xf-xi)/Ns;
w0=3/8;
w1=9/8;
w2=9/8;
w3=3/8;
J=0;
For[i=0,i<=Ns-3,i=i+3,
J=J+dx(f[xi+i dx]w0+f[xi+(i+1) dx]w1+f[xi+(i+2) dx]w2+f[xi+(i+3)dx]w3)];
Return[J];
]


We'll use the integrator above to integrate the following function:


f[x_]:=If[x==0,1,Sin[x]/x];


LagrangeIF[f,0.,1,6]


Using Simpson's method:


Simpson[f_, a_, b_, Ns_] := Module[{Ia, j, dx},
dx = (b-a)/Ns;
Ia = 0.0;
For[j =0, j <= Ns-2, j = j+2,
Ia = Ia + (1/3) dx (f[a + j dx] + 4 f[a + (j+1) dx] + f[a + (j+2) dx]);
];
Return[Ia];
]


Simpson[f,0,1,4]


Problem 6.12


The initial angle will be


\[Theta]i=ArcCos[0.5];


The function to be integrated will be


A=Sin[\[Theta]i/2];
T[\[Phi]_]:=4Sqrt[1/9.8]/Sqrt[1-A^2 (Sin[\[Phi]])^2]


LagrangeIF[T,0,Pi/2]


For the small angle approximation (sin(\[Theta])\[TildeTilde]\[Theta]):


2 Pi Sqrt[1/9.8]


Problem 6.14


Implementing our integrator to calculate the perimeter of an ellipse, we have:


EllipsePL[a_,b_]:=Module[{e,c,Circ},
e=Sqrt[1-b^2/a^2];
c[\[Theta]_]:=4a Sqrt[1-e^2 (Sin[\[Theta]])^2];
Circ=LagrangeIF[c,0,Pi /2];
Return[Circ]
]


And using the Simpson method we have:


EllipsePS[a_,b_]:=Module[{e,c,Circ},
e=Sqrt[1-b^2/a^2];
c[\[Theta]_]:=4a Sqrt[1-e^2 (Sin[\[Theta]])^2];
Circ=Simpson[c,0,Pi /2,4];
Return[Circ]
]


Inum=EllipsePL[2,2];


The difference between the value obtained with the Simpson method and the exact value is:


Abs[Inum-4 Pi]


Our function with a = b = R = 5 will give


EllipsePL[5,5]


With a = 5 and b = 0


Simplify[EllipsePL[5,0]]


With a = 2 and b = 1	


Simplify[EllipsePL[2,1]]


Problem 6.15


Simpson2D[f_,ax_,bx_,Nx_,ay_,by_,Ny_] := Module[{Iy},
Iy[x_] := Module[{Igrandy},
Igrandy[y_] := f[x,y];
Return[Simpson[Igrandy,ay,by,Ny]];
];
Return[Simpson[Iy,ax,bx,Nx]];
]


The equations of the components of the electric field:


Ex[x_,y_,z_]:=Module[{IEx,ex},
IEx[r_,\[Phi]_]:= ((x-r Cos[\[Phi]])r )/((x-r Cos[\[Phi]])^2 +(y-r Sin[\[Phi]])^2+z^2)^(3/2);
ex=Simpson2D[IEx, 0,2,10,0,2 Pi,10];
Return[ex];
]


Ey[x_,y_,z_]:=Module[{IEy,ey},
IEy[r_,\[Phi]_]:= ((y-r Sin[\[Phi]])r )/((x-r Cos[\[Phi]])^2 +(y-r Sin[\[Phi]])^2+z^2)^(3/2);
ey=Simpson2D[IEy, 0,2,10,0,2 Pi,10];
Return[ey];
]


Ez[x_,y_,z_]:=Module[{IEz,ez},
IEz[r_,\[Phi]_]:= (z r )/((x-r Cos[\[Phi]])^2 +(y-r Sin[\[Phi]])^2+z^2)^(3/2);
ez=Simpson2D[IEz, 0,2,10,0,2 Pi,10];
Return[ez];
]


w=Table[{{x,y,z},{Ex[x,y,z],Ey[x,y,z],Ez[x,y,z]}},{x,-5,5,1},{y,-5,5,1},{z,1,5,1}];


e=ListVectorPlot3D[w];
d=Graphics3D[{Purple,Cylinder[{{0,0,0},{0,0,.001}},2]}];
Show[e,d]


The next function will calculate the electric field when x=0 and y=0.


ZField[zi_,zf_,dz_]:=Module[{field,z,i},
field=Table[0.0,{i,1,(zf-zi)+1/dz},{j,1,2}];
z=zi;
i=1;
While[z<=zf,
field[[i,1]]=z;
field[[i,2]]=Ez[0,0,z];
z=z+dz;
i=i+1;
];
Return[field];
]


The theoretical expression for E=E z is:


Ztheoretical[z_]:=2 Pi (1-(z/(Sqrt[4+z^2])))


Plotting the theoretical and the numerical solution together we obtain:


a=ListPlot[ZField[1,10,1]];
b=Plot[Ztheoretical[z],{z,0,10}];
Show[a,b,AxesLabel-> {"z","Ez"}]


Problem 6.17


This function will take a function and approximate it's integral between a lower and an upper boundary for values of N up to a maximum N given by the user.


MonteCarlo[f_,N_,lb_,ub_]:=Module[{dx,n,m,x,I,Idata},
Idata=Table[0.0,{i,1,N}];
m=1;
n=1;
I=0;
While[m<=N,

While[n<= m,
dx=.5/m;
x=RandomReal[{lb,ub}];
I=I+f[x] dx;
n=n+1];

Idata[[m]]=I;
m=m+1;
n=1;
I=0];
Return[Idata];
]


The function whose integral we will approximate is:


f[x_]:=x^5


data=ListPlot[MonteCarlo[f,500,1/2,1],PlotStyle->Red,AxesLabel->{"N","Integral value"}];
tvalue=Plot[21/128,{x,0,500}];
Show[data,tvalue]


The analytical value of the integral is given by the line in the graph, while the dots represent the approximations of the integral for different values of N. We can see that the approximations get closer to the line as N increases.


Problem 6.18


The approximation of the volume of a two dimensional sphere will be given by the next function.


MC2Dsph[R_,N_]:=Module[{n,hit,x,y},
hit=0;
n=1;
While[n<=N,
x=RandomReal[{-R,R}];
y=RandomReal[{-R,R}];
If[Sqrt[x^2+y^2]<=R,
hit=hit+1];
n=n+1;
];
Return[(hit/N) (2R)^2];
]


We will take a sphere of radius 2.


MC2Dsph[2.,1000]


The analytical result for the volume will be:


Pi 2^2.


This function will give the number of "darts" needed to obtain the volume of the two dimensional sphere of radius R up to an error of approximately \[PlusMinus](di/2). The function will take a list with M approximations of the volume and take the difference between the higest and lowest values. It will start with one dart, and will increase the number of darts until the difference between the highest and lowest values is lower or equal to di. 


Hits2D[R_,M_,di_]:=Module[{m,n,Idata},
Idata=Table[0.0,{i,1,M}];
m=1;
n=1;
Idata[[1]]=20;
Idata[[2]]=1;

While[Max[Idata]-Min[Idata]>di,

While[m<=M,
Idata[[m]]=MC2Dsph[R,n];
m=m+1];

n=n+1;
m=1;
];
Return[n];
]


Hits2D[2,5,.06]


Similarly, for a three-dimensional sphere of radius 2:


MC3Dsph[R_,N_]:=Module[{n,hit,x,y,z},
hit=0;
n=1;
While[n<=N,
x=RandomReal[{-R,R}];
y=RandomReal[{-R,R}];
z=RandomReal[{-R,R}];
If[Sqrt[x^2+y^2+z^2]<=R,
hit=hit+1];
n=n+1;
];
Return[(hit/N) (2R)^3];
]


MC3Dsph[2.,500]


(4/3.) Pi 2^3


Hits3D[R_,M_,di_]:=Module[{m,n,Idata},
Idata=Table[0.0,{i,1,M}];
m=1;
n=1;
Idata[[1]]=20;
Idata[[2]]=1;

While[Max[Idata]-Min[Idata]>di,

While[m<=M,
Idata[[m]]=MC3Dsph[R,n];
m=m+1];

n=n+1;
m=1;
];
Return[n];
]


Hits3D[2,3,0.01]


Finally, for a four-dimensional sphere of radius 2:


MC4Dsph[R_,N_]:=Module[{n,hit,x,y,z,u},
hit=0;
n=1;
While[n<=N,
x=RandomReal[{-R,R}];
y=RandomReal[{-R,R}];
z=RandomReal[{-R,R}];
u=RandomReal[{-R,R}];
If[Sqrt[x^2+y^2+z^2+u^2]<=R,
hit=hit+1];
n=n+1;
];
Return[(hit/N) (2R)^4];
]


MC4Dsph[2.,700]


(1/2.) Pi^2 2^4 


Hits4D[R_,M_,di_]:=Module[{m,n,Idata},
Idata=Table[0.0,{i,1,M}];
m=1;
n=1;
Idata[[1]]=20;
Idata[[2]]=1;

While[Max[Idata]-Min[Idata]>di,

While[m<=M,
Idata[[m]]=MC4Dsph[R,n];
m=m+1];

n=n+1;
m=1;
];
Return[n];
]


Hits4D[2,3,0.09]
