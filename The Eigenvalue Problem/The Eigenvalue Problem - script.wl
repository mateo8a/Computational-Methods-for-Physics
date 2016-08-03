(* ::Package:: *)

Lab 9


Mateo Ochoa


Problem 10.11


data=Import["http://people.reed.edu/~jfrankli/Book/downloads/files/ch10.dat"];


FFT[invec_] := Module[{k,elist,olist,netlist,eft,oft,n,omega},
If[Length[invec] ==  1, 
Return[invec];
];
n = Length[invec];
omega = 1;
elist = Table[invec[[j]],{j,1,n,2}];
olist = Table[invec[[j]],{j,2,n,2}];
eft = FFT[elist];
oft = FFT[olist];
netlist = Table[0,{j,0,n-1}];
For[k = 0, k <= n/2 - 1, k = k+1,
netlist[[k+1]] = eft[[k+1]] + omega oft[[k+1]];
netlist[[k + n/2+1]] = eft[[k+1]] - omega oft[[k+1]];
omega = omega Exp[2 Pi I/n];
];
Return[netlist];
]


psubj=Table[data[[i,2]],{i,1,Length[data]}];


n=Length[psubj];
dt=data[[2,1]]-data[[1,1]];
df = 1/(Length[psubj] dt);
FTTtransform=FFT[psubj];
FTTtransform=RotateRight[FTTtransform,n/2];
pspec=Table[{k df,(Abs[FTTtransform[[(k+n/2)+1]]])^2},{k,-n/2,n/2 -1}];


ListLinePlot[pspec,PlotRange->{{-3,3},All},AxesLabel->{"\!\(\*SubscriptBox[\(f\), \(k\)]\)","|P(\!\(\*SubscriptBox[\(f\), \(k\)]\))\!\(\*SuperscriptBox[\(|\), \(2\)]\)"}]


Fk[data_]:=Module[{i,list,listf,k,kk,maxv},
list=Table[0,{i,1,Length[data]}];
k=1;
maxv=Max[(Abs[data])^2];
For[i=1,i<=Length[data],i=i+1,
If[(Abs[data[[i]]])^2==maxv,list[[k]]={i-1,data[[i]]};k=k+1]
];
kk=k-1;
listf=Table[0,{i,1,kk}];
For[i=1,i<=kk,i=i+1,
listf[[i]]=list[[i]]];
Return[listf]]


fk=Fk[FTTtransform];


freqHz[data_]:=Module[{i,list,listf,k,kk},
list=Table[0,{i,1,Length[data]}];
k=1;
For[i=1,i<=Length[data],i=i+1,
list[[k]]={(-n df/2)+data[[i,1]] df,data[[i,2]]};k=k+1];
kk=k-1;
listf=Table[0,{i,1,kk}];
For[i=1,i<=kk,i=i+1,
listf[[i]]=list[[i]]];
Return[listf]]


The frequency of oscillation of the mass is:


fnc=Chop[freqHz[fk]];
MatrixForm[fnc]


UDBacksolve[U_,y_] := Module[{xout,i,j,NN},
NN = Length[U];
xout = Table[0,{ii,1,NN}];
xout[[NN]] = y[[NN]]/U[[NN,NN]];
For[i = NN-1, i >= 1, i = i-1,
xout[[i]] = N[(y[[i]]- Sum[U[[i,jj]]xout[[jj]] ,{jj,i+1,NN}])/U[[i,i]]];
];
Return[xout];
]


The functions that we will use to fit the data have the form:


x1[t_,b_]:=Exp[-b t] (Cos[w  t])


x2[t_,b_]:=Exp[-b t]((b/w)Sin[w  t])


w=2 Pi fnc[[2,1]];


tdata=Table[data[[i,1]],{i,1,Length[data]}];


This function will generate a matrix with a set of functions evaluated at the values of tdata.


FuncMatrix[f1_,f2_,b_,tdata_]:=Module[{matrix,bb,tt,n,tinput},
n=Length[tdata];
matrix=Table[0.0,{i,1,n},{j,1,2*b}];
For[bb=1,bb<=b,bb=bb+1,
For[tt=1,tt<=n,tt=tt+1,
tinput=tdata[[tt]];
matrix[[tt,bb]]=Chop[f1[tinput,bb]];
];
];
For[bb=b+1,bb<=2b,bb=bb+1,
For[tt=1,tt<=n,tt=tt+1,
tinput=tdata[[tt]];
matrix[[tt,bb]]=Chop[f2[tinput,bb-b]];
];
];
Return[matrix];
]


Let our matrix A be:


A=FuncMatrix[x1,x2,4,tdata];


with evaluated functions that have b with values up to 4.


{Qq,Rr}=QRDecomposition[A];


xdata=Table[data[[i,2]],{i,1,Length[data]}];


y=Qq.xdata;


The coefficients \[Alpha] will be:


coeff=UDBacksolve[Rr,y];


r^2 will be:


r=xdata-A.coeff;
rsquared=r.r;


The variance of the data:


averagex=Sum[xdata[[i]],{i,1,Length[data]}]/Length[data];
variance=Sum[(xdata[[i]]-averagex)^2,{i,1,Length[xdata]}]/Length[data];


Finally, R^2:


Rsquared=1-rsquared/variance


Problem 10.13


The modified Gram-Schmidt procedure is:


NZERO = 10^(-10);


GSmod[A_] := Module[{Q,R,i,j,At,qi},
At = Transpose[A];
R = Table[0,{ii,1,Length[A]},{jj,1,Length[A]}];
Q = {};
For[j = 1, j <= Length[A], j = j+1,
qi =N[ At[[j]] ];
For[i = 1, i <= j - 1, i = i+1,
R[[i,j]] = Q[[i]].qi;
qi = qi - R[[i,j]] Q[[i]];
];
R[[j,j]] = Sqrt[qi.qi];
If[Abs[R[[j,j]]] <= NZERO, Print["Matrix A may be singular"]; Return[{Null,Null}];];
Q = Join[Q,N[{qi/R[[j,j]]}]];
];
Return[{Transpose[Q],R}]
]


The matrix to be used is:


A={{1,2,3,4},{2,5,6,7},{3,6,8,9},{4,7,9,10}};


SIsolver[matrix_,M_]:=Module[{Qo,j,A,Z,Qj,Rj,QR,Aj},
A=matrix;
Qo=Table[If[i==j,1,0],{i,1,Length[A]},{j,1,Length[A]}];
Qj=Qo;
For[j=1,j<=M,j=j+1,
Z=A.Qj;
QR=GSmod[Z];
Qj=QR[[1]];
Rj=QR[[2]];
Aj=Transpose[Qj].A.Qj;
];
Return[{Qj,Aj}];
]


l=SIsolver[A,100];


MatrixForm[l[[1]]]


MatrixForm[Chop[l[[2]]]]


QRsolver[matrix_,M_]:=Module[{j,Qj,Rj,QR,Aj,Vo,Vj},
Aj=matrix;
Vo=Table[If[i==j,1,0],{i,1,Length[Aj]},{j,1,Length[Aj]}];
Vj=Vo;
For[j=1,j<=M,j=j+1,
QR=GSmod[Aj];
Qj=QR[[1]];
Rj=QR[[2]];
Aj=Rj.Qj;
Vj=Vj.Qj;
];
Return[{Vj,Aj}];
]


l=QRsolver[A,100];


MatrixForm[l[[1]]]


MatrixForm[Chop[l[[2]]]]


The following function finds the number of iterations M necessary for all the off-diagonal entries of the Aj matrix to have a value less than an specified eps. One can specify the function to be used to find the eigenvalues and eigenvectors.


Mfinder[matrix_,eps_,method_]:=Module[{A,M,QA,Aj,entry,l,i,j,elements,n},
A=matrix;
elements=Table[0.0,{i,1,Length[A]*Length[A]}];
l=Length[A];
M=1;
entry=eps+1;
n=0;
While[entry>eps,
QA=method[A,M];
Aj=QA[[2]];
For[i=1,i<=l,i=i+1,
For[j=1,j<=l,j=j+1,
If[j==i,
n=n+1;
elements[[n]]=0,
n=n+1;
elements[[n]]=Aj[[i,j]];
];
];
];
entry=Abs[Sum[elements[[i]],{i,Length[elements]}]/(Length[elements]-Length[Aj])];
M=M+1;
n=0;
];
Return[M];
]


For the SI iteration:


Mfinder[A,10^(-13),SIsolver]


For the QR iteration:


Mfinder[A,10^(-13),QRsolver]


Problem 10.15


matrixD[Nq_, qinf_] := Module[{dq,matrix,i},
dq=qinf/(Nq+1);

matrix = Table[0.0,{i,1,Nq},{j,1,Nq}];
matrix[[1,1 ]] += -2.0/(-dq^2) -2/dq;
matrix[[1,2 ]] += 1.0/(-dq^2);
matrix[[Nq,Nq]] += -2.0/(-dq^2) -2/(dq Nq);
matrix[[Nq,Nq-1]] += 1.0/(-dq^2);

For[i=2,i<=Nq-1,i=i+1,
matrix[[i,i-1]] += 1.0/(-dq^2);
matrix[[i,i]] += -2.0/(-dq^2) -2/(i dq);
matrix[[i,i+1]] += 1.0/(-dq^2);
];
matrix=matrix;
Return[matrix];
]


Nq=10000;
qinf=100;
DD=matrixD[Nq,qinf];


This function will obtain the first n negative eigenvalues of the desired matrix.


Evals[matrix_,n_]:=Module[{eigenvals,evalsmat,evalsmat2,i,k},
evalsmat2=Table[0.0,{i,1,n}];
eigenvals=Eigenvalues[matrix];
evalsmat=Table[0.0,{i,1,Length[eigenvals]}];
k=1;

For[i=1,i<=Length[eigenvals],i=i+1,
If[eigenvals[[i]]<0,
evalsmat[[k]]=eigenvals[[i]];
k=k+1;
];
];
For[i=1,i<= n,i=i+1,
evalsmat2[[i]]=evalsmat[[i]]];
Return[evalsmat2];
]


evals=Evals[DD,6]


The first six negative bound-state energies of hydrogen are:


MatrixForm[evals]


The functions that were used to approximate the data have the form:


f[n_,exp_]:=1/n^exp


We know that n will go from 1 to 6.


ndata={1,2,3,4,5,6};


FuncMatrix[f1_,exp_,ndata_]:=Module[{matrix,bb,nn,NN,ninput},
NN=Length[ndata];
matrix=Table[0.0,{i,1,NN},{j,1,exp}];
For[bb=1,bb<=exp,bb=bb+1,
For[nn=1,nn<=NN,nn=nn+1,
ninput=ndata[[nn]];
matrix[[nn,bb]]=Chop[f1[ninput,bb]];
];
];
Return[matrix];
]


A=FuncMatrix[f,3,ndata];


{Qq,Rr}=QRDecomposition[A];


Edata=evals;


y=Qq.Edata;


coeff=UDBacksolve[Rr,y];


r=Edata-A.coeff;
rsquared=r.r;


averageE=Sum[Edata[[i]],{i,1,Length[Edata]}]/Length[Edata];
variance=Sum[(Edata[[i]]-averageE)^2,{i,1,Length[Edata]}]/Length[Edata];


Rsquared=1-rsquared/variance


Problem 10.16


data=Import["http://people.reed.edu/~jfrankli/Book/downloads/files/ch10.mystery.mat","csv"];


eigenv=Eigensystem[data];


Taking n=4:


n=4;
Vs=Take[eigenv[[2]],{1,n}];
Es=DiagonalMatrix[Take[eigenv[[1]],{1,n}]];
matrix=Transpose[Vs].Es.Vs;
MatrixPlot[matrix]


Taking n=100:


n=100;
Vs=Take[eigenv[[2]],{1,n}];
Es=DiagonalMatrix[Take[eigenv[[1]],{1,n}]];
matrix=Transpose[Vs].Es.Vs;
MatrixPlot[matrix]
