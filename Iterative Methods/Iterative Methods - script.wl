(* ::Package:: *)

Lab 10


Mateo Ochoa


Problem 11.12


This function will generate three matrices, one with the diagonal entries, and two strictly upper and lower diagonal portions of the input matrix.


DLUdecomp[matrix_]:=Module[{D,U,L,nn,i,j},
nn=Length[matrix];
D=Table[0.0,{i,1,nn},{i,1,nn}];
U=Table[0.0,{i,1,nn},{i,1,nn}];
L=Table[0.0,{i,1,nn},{i,1,nn}];
For[i=1,i<=nn,i=i+1,
For[j=1,j<=nn,j=j+1,
If[i==j,D[[i,j]]=matrix[[i,j]]];
];
];
For[j=2,j<=nn,j=j+1,
For[i=1,i<j,i=i+1,
U[[i,j]]=matrix[[i,j]];
];
];
For[j=1,j<nn,j=j+1,
For[i=nn,i>j,i=i-1,
L[[i,j]]=matrix[[i,j]];
];
];
Return[{D,L,U}];
]


Chapter 4 functions


DXmat[Nx_, Ny_, dx_, dy_, left_, right_] := Module[{g,retmat,bside, index, indey,i,j,k},
g[n_,m_] := (m - 1)*Nx + n;
retmat=SparseArray[{},{Nx Ny,Nx Ny},0];
bside=SparseArray[{},{Nx Ny},0];
For[indey = 1, indey <= Ny, indey = indey+1,
j = g[1,indey];
k = g[2,indey];
retmat= retmat+SparseArray[{{j,j}->-2.0/dx^2,{j,k}->1.0/dx^2},{Nx Ny,Nx Ny},0];
bside=bside+SparseArray[{j->-left[indey dy]/dx^2},{Nx Ny},0];
For[index = 2, index <= Nx-1, index = index+1,
i = g[index-1,indey];
j = g[index,indey];
k = g[index + 1,indey];
retmat=retmat+SparseArray[{{j,i}->1.0/dx^2,{j,j}->-2.0/dx^2,{j,k}->1.0/dx^2},{Nx Ny,Nx Ny},0];
];
i = g[Nx-1,indey];
j = g[Nx,indey];
retmat= retmat+SparseArray[{{j,i}->1.0/dx^2,{j,j}->-2.0/dx^2},{Nx Ny,Nx Ny},0];
bside=bside+SparseArray[{j->-right[indey dy]/dx^2},{Nx Ny},0];
];
Return[{retmat, bside}];
]


DYmat[Nx_, Ny_, dx_, dy_, top_, bottom_] := Module[{g,retmat,bside, index, indey,i,j,k},
g[n_,m_] := (m - 1)*Nx + n;
retmat=SparseArray[{},{Nx Ny,Nx Ny},0];
bside=SparseArray[{},{Nx Ny},0];
For[index = 1, index <= Nx, index = index+1,
j = g[index,1];
k = g[index, 2];
retmat= retmat+SparseArray[{{j,j}->-2.0/dy^2,{j,k}->1.0/dy^2},{Nx Ny,Nx Ny},0];
bside=bside+SparseArray[{j->-bottom[indey dx]/dy^2},{Nx Ny},0];
For[indey = 2, indey <= Ny-1, indey = indey+1,
i = g[index,indey-1];
j = g[index,indey];
k = g[index ,indey+1];
retmat=retmat+SparseArray[{{j,i}->1.0/dy^2,{j,j}->-2.0/dy^2,{j,k}->1.0/dy^2},{Nx Ny,Nx Ny},0];
];
i = g[index, Ny-1];
j = g[index, Ny];
retmat= retmat+SparseArray[{{j,i}->1.0/dy^2,{j,j}->-2.0/dy^2},{Nx Ny,Nx Ny},0];
bside=bside+SparseArray[{j->-top[indey dx]/dy^2},{Nx Ny},0];
];
Return[{retmat, bside}];
]


MakeSource[ Nx_, Ny_, dx_, dy_, source_] := Module[{svec, g,j, index, indey},
g[n_,m_] := (m - 1)*Nx + n;
svec = SparseArray[{},{Nx Ny},0];
For[index = 1, index <= Nx, index=index+1,
For[indey = 1, indey <= Ny, indey = indey + 1,
j = g[index,indey];
svec = svec+SparseArray[{j->source[ index dx, indey dy]},{Nx Ny},0];
];
];
Return[svec];
]


Part (a)


RHS[x_,y_] :=Module[{ret},
ret = 0.0;
If[Round[x/dx]==15 && Round[y/dy]==15,ret=1];
Return[ret];
]


top[x_] :=0.;
bottom[x_] :=0.;
left[y_] := 0.;
right[y_] :=0.;


Nx = 29;
Ny = 29;
dx = 1/(Nx + 1);
dy = 1/(Ny+1);


{DX, bx} =DXmat[Nx, Ny,dx,dy,  left, right] ;
{DY, by} = DYmat[Nx, Ny,  dx,dy,top, bottom];
bs = MakeSource[ Nx, Ny, dx, dy,RHS] ;
b = bx + by + bs;
A = DX + DY;


w=1.9;
n=300.;


The SOR method will be implemented as follows.


SOR[A_,b_,n_,w_]:=Module[{DD,LL,UU,l,u,i,x,rhs,lhs,xo},
SeedRandom[1];
xo=Table[RandomReal[{-1,1}],{i,1,Nx Ny}];
x=xo;
{DD,l,u}=DLUdecomp[A];
LL=-l;
UU=-u;
For[i=1,i<=n,i=i+1,
rhs=((1-w)*DD+w* UU).x+w b;
lhs=DD-w LL;
x=LinearSolve[lhs,rhs];
];
Return[x];
]


The solution after 300 iterations is:


x=SOR[A,b,n,w];


Where the maximum difference between the ith elements of Ax-b is:


Max[A.x-b]


Plotting this, we have:


g[n_,m_] := (m - 1)*Nx + n;
Datatable= Table[{ii dx, jj dy, x[[ g[ii,jj] ]]},{ii,1, Nx}, {jj,1,Ny}];
Datatable = Partition[Flatten[Datatable],3];
Ga = ListContourPlot[Datatable,Contours-> 20,ColorFunction-> "GrayTones"]


Part (b)


MaxRes[A_,b_,w_,mi_,mf_]:=Module[{maxres,i,x,max},
maxres=Table[0.0,{i,1,mf-mi+1}];
For[i=mi,i<=mf,i=i+1,
x=SOR[A,b,i,w];
max=Max[A.x-b];
maxres[[i-mi+1]]=max;
];
Return[maxres];
]


w=0.5;
maxr=MaxRes[A,b,w,10,300];
maxr=Abs[maxr];
ListPlot[maxr]


w=1;
maxr=MaxRes[A,b,w,10,300];
maxr=Abs[maxr];
ListPlot[maxr]


w=1.5;
maxr=MaxRes[A,b,w,10,300];
maxr=Abs[maxr];
ListPlot[maxr]


Problem 11.14


matrixD[Nq_, qinf_,Eo_] := Module[{dq,matrix,i},
dq=qinf/(Nq+1);

matrix = Table[0.0,{i,1,Nq},{j,1,Nq}];
matrix[[1,1 ]] += -2.0/(-dq^2) +Eo-2/dq;
matrix[[1,2 ]] += 1.0/(-dq^2);
matrix[[Nq,Nq]] += -2.0/(-dq^2) +Eo-2/(dq Nq);
matrix[[Nq,Nq-1]] += 1.0/(-dq^2);

For[i=2,i<=Nq-1,i=i+1,
matrix[[i,i-1]] += 1.0/(-dq^2);
matrix[[i,i]] += -2.0/(-dq^2) +Eo-2/(i dq);
matrix[[i,i+1]] += 1.0/(-dq^2);
];
matrix=matrix;
Return[matrix];
]


The set-up for the problem will be:


Nq=1000;
qinf=40;
Eo=1;
DD=matrixD[Nq,qinf,Eo];


Lanczos[inA_, inb_,m_] := Module[{retQ,retT,alpha,beta,n,w,v,index},
n = Length[inA];
retQ = Table[Table[0.0,{i,1,n}],{j,1,m+2}];
alpha = Table[0.0,{i,1,m+1}];
beta = Table[0.0,{j,1,m+1}];
retQ[[2]] = inb/Sqrt[inb.inb];
For[index = 2, index <=  m+1, index=index+1,
w = inA.retQ[[index]];
alpha[[index]] = retQ[[index]].w;
v = w - beta[[index-1]] retQ[[index-1]] - alpha[[index]] retQ[[index]];
retQ[[index+1]] = v/Sqrt[v.v];
beta[[index]] = retQ[[index+1]].v;
If[Abs[beta[[index]]] <= 10^(-10), Print["Warning - try another vector b"];];
];
retQ = retQ[[2;;m+1]];
beta = beta[[2;;m]];
alpha = alpha[[2;;m+1]];
retT = DiagonalMatrix[alpha] + DiagonalMatrix[beta,1] + DiagonalMatrix[beta,-1];
 Return[{Transpose[retQ],retT}]; 
]


The initial vector will be a random vector with norm 1.


b=Table[RandomReal[{-50,50}],{i,1,Length[Transpose[DD]]}];
bb=Sqrt[b.b];
b=b/bb;
m=300;


Dmatrix=Lanczos[DD,b,m];


In order to find the bound states of hydrogen, the value of Eo is substracted from the eigenvalues of the matrix T given by the Lanczos method in the next function.


Evals[matrix_,Eo_]:=Module[{evals,evtable,evtable2,i,k},
evals=Eigenvalues[Dmatrix[[2]]]-Eo;
evtable=Table[0.0,{i,1,Length[evals]}];
k=1;
For[i=1,i<=Length[evals],i=i+1,
If[evals[[i]]<0,evtable[[k]]=evals[[i]];k=k+1];
];
evtable2=Table[0.0,{i,1,k-1}];
For[i=1,i<=k-1,i=i+1,
evtable2[[i]]=evtable[[i]]];
Return[evtable2];
]


The approximate eigenvalues of the matrix A, which correspond to the energy spectrum of hydrogen will be:


Evals[Dmatrix[[2]],Eo]


Problem 11.15


matrixD[Nq_, qinf_,Eo_] := Module[{dq,matrix,i},
dq=qinf/(Nq+1);

matrix = Table[0.0,{i,1,Nq},{j,1,Nq}];
matrix[[1,1 ]] += -2.0/(-dq^2) +Eo-2/dq;
matrix[[1,2 ]] += 1.0/(-dq^2);
matrix[[Nq,Nq]] += -2.0/(-dq^2) +Eo-2/(dq Nq);
matrix[[Nq,Nq-1]] += 1.0/(-dq^2);

For[i=2,i<=Nq-1,i=i+1,
matrix[[i,i-1]] += 1.0/(-dq^2);
matrix[[i,i]] += -2.0/(-dq^2) +Eo-2/(i dq);
matrix[[i,i+1]] += 1.0/(-dq^2);
];
matrix=matrix;
Return[matrix];
]


The set-up for the problem will be:


Nq=1000;
qinf=40;
Eo=1;
A=matrixD[Nq,qinf,Eo];


The, we create a matrix B \[Congruent] (A - \[Mu] (I))^-1:


mu=Eo-1+1/100;
B=Inverse[A-mu*IdentityMatrix[Length[A]]];


Then, we obtain the eigenvalues of the matrix B. 


b=Table[RandomReal[{-50,50}],{i,1,Length[Transpose[DD]]}];
bb=Sqrt[b.b];
b=b/bb;
m=4;
Dmatrix=Lanczos[B,b,m];
Ev=Evals[Dmatrix[[2]],Eo]


Comparing to the values of the grond state from the previous lab, 


we obtain the error of the newly found ground state.


error=100*Abs[(Ev[[3]]+0.9975173655152572`)]/0.9975173655152572`


The new approximation is very accurate.


Solving again for a new value of \[Mu] that will help us obtain the third eigenvalue of the hydrogen spectrum, we have:


mu=Eo-(1/9)+1/100;
B=Inverse[A-mu*IdentityMatrix[Length[A]]];


b=Table[RandomReal[{-50,50}],{i,1,Length[Transpose[DD]]}];
bb=Sqrt[b.b];
b=b/bb;
m=4;
Dmatrix=Lanczos[B,b,m];
Ev=Evals[Dmatrix[[2]],Eo]
