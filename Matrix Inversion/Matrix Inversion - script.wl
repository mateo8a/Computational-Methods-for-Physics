(* ::Package:: *)

Lab 8


Mateo Ochoa


9.10


LUdecomp[matrixA_]:=Module[{n,LL,UU,k,j,Lk,Linvk,det,jj},
If[Length[matrixA]==1,Print["Error. Singular matrix."],
n=Length[matrixA];
LL=Table[If[i==j,1,0],{i,1,n},{j,1,n}];
UU=matrixA;
For[k=1,k<= n-1,k=k+1,
Lk=Table[If[i==j,1,0],{i,1,n},{j,1,n}];
Linvk=Table[If[i==j,1,0],{i,1,n},{j,1,n}];
For[j=k+1,j<=n,j=j+1,
Lk[[j,k]]=-UU[[j,k]]/UU[[k,k]];
Linvk[[j,k]]=-Lk[[j,k]];
];
UU=Lk.UU;
LL=LL.Linvk;
];
det=1;
For[jj=1,jj<=n,jj=jj+1,
det=det*UU[[jj,jj]]*LL[[jj,jj]]
];
Return[det];
];
]


This function uses LU decomposition to construct the determinant of the input matrix.


A={{1,2,3},{4,5,6},{7,8,10}};


LUdecomp[A]


We defined the matrix A as above, and obtained its determinant.


Problem 9.11


LUdecomp[matrixA_]:=Module[{n,LL,UU,k,j,Lk,Linvk,det,jj},
If[Length[matrixA]==1,Print["Error. Singular matrix."],
n=Length[matrixA];
LL=Table[If[i==j,1,0],{i,1,n},{j,1,n}];
UU=matrixA;
For[k=1,k<= n-1,k=k+1,
Lk=Table[If[i==j,1,0],{i,1,n},{j,1,n}];
Linvk=Table[If[i==j,1,0],{i,1,n},{j,1,n}];
For[j=k+1,j<=n,j=j+1,
Lk[[j,k]]=-UU[[j,k]]/UU[[k,k]];
Linvk[[j,k]]=-Lk[[j,k]];
];
UU=Lk.UU;
LL=LL.Linvk;
];
det=1;
For[jj=1,jj<=n,jj=jj+1,
det=det*UU[[jj,jj]]*LL[[jj,jj]]
];
Return[{UU,LL}];
];
]


This function uses LU decomposition to construct the upper and lower triangular matrices corresponding to the input matrix.


UDBacksolve[U_,y_] := Module[{xout,i,j,NN},
NN = Length[U];
xout = Table[0,{ii,1,NN}];
xout[[NN]] = y[[NN]]/U[[NN,NN]];
For[i = NN-1, i >= 1, i = i-1,
xout[[i]] = N[(y[[i]]- Sum[U[[i,jj]]xout[[jj]] ,{jj,i+1,NN}])/U[[i,i]]];
];
Return[xout];
]


This function solves for x in Ux = y, where U is an upper triangular matrix.


LDBacksolve[L_,b_] := Module[{yout,i,j,NN},
NN=Length[L];
yout=Table[0,{ii,1,NN}];
yout[[1]]=b[[1]]/L[[1,1]];
For[i=2,i<=NN,i=i+1,
yout[[i]]=N[(b[[i]]-Sum[L[[i,jj]] yout[[jj]],{jj,1,i-1}])/L[[i,i]]];
];
Return[yout];
]


This function solves for y in Ly = b, where L is a lower triangular matrix.


LUSolver[A_,b_]:=Module[{Umatrix,Lmatrix,matrices,yvector,xvector},
matrices=LUdecomp[A];
Umatrix=matrices[[1]];
Lmatrix=matrices[[2]];
yvector=LDBacksolve[Lmatrix,b];
xvector=UDBacksolve[Umatrix,yvector];
Return[xvector];
]


This function solves for x in Ax = b by converting A into upper and lower triangular matrices,and using the functions UDBacksolve and LDBacksolve.


residual[A_,b_,x_]:=Module[{residual,max},
residual=A.x - b;
max=Max[Abs[residual]];
Return[max];
]


This function calculates the residual r = Ax - b.


Ninversions[N_,dimension_,lowerb_,upperb_]:=Module[{i,n,xvector,maxresidual,xlist,A,b,AA,bb,ret},
xlist=Table[0,{i,1,N}];
AA=Table[0,{i,1,N}];
bb=Table[0,{i,1,N}];
ret=Table[0,{i,1,N}];
n=dimension;
For[i=1,i<=N,i=i+1,
A=Table[RandomReal[{lowerb,upperb}],{i,1,n},{j,1,n}];
b=Table[RandomReal[{lowerb,upperb}],{i,1,n}];
xvector=LUSolver[A,b];
xlist[[i]]=xvector;
AA[[i]]=A;
bb[[i]]=b;
];
ret={AA,bb,xlist};
Return[ret];
]


This function takes the number of inversions N to be done to the matrix whose entries are randomly generated (the values of the entries will be between lowerb and upperb). The dimension input indicates the dimesion of the matrix A and the vector b.


resNinv[N_,dimension_,lowerb_,upperb_]:=Module[{val,i,res,A,b,x},
val=Ninversions[N,dimension,lowerb,upperb];
res=Table[0,{i,1,N}];
For[i=1,i<=N,i=i+1,
A=val[[1,i]];
b=val[[2,i]];
x=val[[3,i]];
res[[i]]=residual[A,b,x];
];
Return[res];
]


This function outputs the residual for N inversions.


Our variables will be:


n=100;
dim=100;
lowb=-1;
upb=1;


residuals=resNinv[n,dim,lowb,upb];
average=Sum[residuals[[i]],{i,1,Length[residuals]}]/Length[residuals];
avplot=Plot[average,{x,0,Length[residuals]+1},PlotStyle->Gray];
resplot=ListPlot[residuals];
Show[avplot,resplot]


resplot=ListPlot[residuals,PlotRange->All]


Problem 9.12


n=100;
A=Table[RandomReal[{-1,1}],{i,1,n},{j,1,n}];


Minverse[Matrix_]:=Module[{i,Matinv,M},
xvector=Table[0.0,{i,1,Length[Matrix]}];
e=Table[If[i==j,1,0],{i,1,Length[Matrix]},{j,1,Length[Matrix]}];
For[i=1,i<=Length[Matrix],i=i+1,
xvector[[i]]=LUSolver[A,e[[i]]];
];
Matinv=Table[xvector[[i]],{i,1,Length[Matrix]}];
M=Transpose[Matinv];
Return[M];
]


Ainv=Minverse[A];


Since n=100 is to show, let n=10 to show that Minverse actually gves the inverse.


n=10;
A=Table[RandomReal[{-1,1}],{i,1,n},{j,1,n}];


Ainv=Minverse[A];
Chop[A.Ainv];
MatrixForm[%]


We can verify that the inverse obtained for A is correct because when the matrix is multiplied by its inverse the identity matrix is returned.


Problem 9.18


Chapter 4 functions


DXmat[Nx_, Ny_, dx_, dy_, left_, right_] := Module[{g,retmat,bside, index, indey,i,j,k},
g[n_,m_] := (m - 1)*Nx + n;

retmat = Table[0.0,{i,1,Nx Ny},{j,1,Nx Ny}];
bside = Table[0.0,{i,1,Nx Ny}];

For[indey = 1, indey <= Ny, indey = indey+1,
j = g[1,indey];
k = g[2,indey];

retmat[[j,j ]] += -2.0/dx^2;
retmat[[j,k ]] += 1.0/dx^2;
bside[[ j ]] += -left[indey dy]/dx^2;

For[index = 2, index <= Nx-1, index = index+1,
i = g[index-1,indey];
j = g[index,indey];
k = g[index + 1,indey];
retmat[[j,i]] += 1.0/dx^2;
retmat[[j,j]] += -2.0/dx^2;
retmat[[j,k]] += 1.0/dx^2;
];
i = g[Nx-1,indey];
j = g[Nx,indey];
retmat[[j,i ]] += 1.0/dx^2;
retmat[[ j,j]] += -2.0/dx^2;
bside[[ j]] += -right[indey dy]/dx^2;
];
Return[{retmat, bside}];
]


DYmat[Nx_, Ny_, dx_, dy_, top_, bottom_] := Module[{g,retmat,bside, index, indey,i,j,k},
g[n_,m_] := (m - 1)*Nx + n;

retmat = Table[0.0,{i,1,Nx Ny},{j,1,Nx Ny}];
bside = Table[0.0,{i,1,Nx Ny}];

For[index = 1, index <= Nx, index = index+1,
j = g[index,1];
k = g[index, 2];

retmat[[j,j ]] += -2.0/dy^2;
retmat[[j,k ]] += 1.0/dy^2;
bside[[ j ]] += -bottom[index dx]/dy^2;

For[indey = 2, indey <= Ny-1, indey = indey+1,
i = g[index,indey-1];
j = g[index,indey];
k = g[index ,indey+1];
retmat[[j,i]] += 1.0/dy^2;
retmat[[j,j]] += -2.0/dy^2;
retmat[[j,k]] += 1.0/dy^2;
];
i = g[index, Ny-1];
j = g[index, Ny];
retmat[[j,i ]] += 1.0/dy^2;
retmat[[ j,j]] += -2.0/dy^2;
bside[[ j]] += -top[index dx]/dy^2;
];
Return[{retmat, bside}];
]



MakeSource[ Nx_, Ny_, dx_, dy_, source_] := Module[{svec, g,j, index, indey},
g[n_,m_] := (m - 1)*Nx + n;

svec = Table[0.0,{i,1,Nx Ny}];
For[index = 1, index <= Nx, index=index+1,
For[indey = 1, indey <= Ny, indey = indey + 1,
j = g[index,indey];
svec[[j ]] = source[ index dx, indey dy];
];
];
Return[svec];
]


top[x_] := 0;
bottom[x_] := 0;
left[y_] := 0;
right[y_] := 0;


Grid definition:


Nx =20;
Ny = 20;
dx = 1/(Nx+ 1);
dy = 1/(Ny+1);


RHS[x_,y_] := Module[{retval, rho},
rho = 1;
retval = 0;
If[y == Round[Ny/4] dy && (x <=(Nx -  Round[Nx/5]+1) dx && x >= (Round[Nx/5] dx)) , 
retval = -rho;
];
If[y == (Ny - Round[Ny/4]+1) dy && (x <=(Nx -  Round[Nx/5]+1) dx && x >= (Round[Nx/5] dx)) , 
retval = rho;
];
Return[-retval];
]


{DX, bx} = DXmat[Nx, Ny, dx, dy, left, right] ;
{DY, by} = DYmat[Nx, Ny, dx, dy, top, bottom];
bs = MakeSource[ Nx, Ny, dx, dy,RHS] ;
rhsvec = bx + by + bs;
D2mat = DX + DY;


GradV[Vxy_,Nx_, Ny_, dx_,dy_,sign_] := Module[{g, gradout, i, j,gv,ij,hp1,hm1,vp1,vm1},
g[n_,m_] := (m - 1)*Nx + n;

gradout = {};
For[i = 2, i <= Nx - 1, i=i+1,
For[j = 2, j<= Ny - 1, j = j + 1,
hp1 = g[i+1,j];
hm1 = g[i-1,j];
vp1 = g[i,j+1];
vm1 = g[i,j-1];
gv ={ (Vxy[[hp1]] - Vxy[[hm1]])/(2 dx),(Vxy[[vp1]] - Vxy[[vm1]])/(2 dy)};
gradout = Join[gradout,{{{i dx, j dy},sign gv}}];
];
];
Return[gradout];
]


Solve using LU decomposition


fxy = LUSolver[D2mat, rhsvec];


g[n_,m_] := (m - 1)*Nx + n
Datatable= Table[{ii dx, jj dy, fxy[[ g[ii,jj] ]]},{ii,1, Nx}, {jj,1,Ny}];
Datatable = Partition[Flatten[Datatable],3];
Ga = ListContourPlot[Datatable,Contours-> 20,ColorFunction-> "GrayTones"]	


Efield = N[GradV[fxy,Nx, Ny,dx,dy,1]];
GE = ListVectorPlot[Efield, VectorStyle->GrayLevel[0.5]]
