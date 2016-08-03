(* ::Package:: *)

Lab 4


Mateo Ochoa


Problem 4.11


The functions that will be used to generate the matrices Dx and Dy will be:


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



The point source at n = m = 25 will be given by:


RHS[x_,y_] :=Module[{ret},
ret = 0.0;
If[Round[x/dx]==25 && Round[y/dy]==25,ret=1];
Return[ret];
]


The source vector will be given by:


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


The boundary values will be:


top[x_] :=0.


bottom[x_] :=0.


left[y_] := 0.


right[y_] :=0.


The matrix will be 49 x 49:


Nx = 49;
Ny = 49;
dx = 1/(Nx + 1);
dy = 1/(Ny+1);


Solving for the derivative matrices DX and DY, the boundary vectors bx and by, the source vector bs, and putting them together (with the elements in their respective sides of the equation) we have:


{DX, bx} =DXmat[Nx, Ny,dx,dy,  left, right] ;
{DY, by} = DYmat[Nx, Ny,  dx,dy,top, bottom];
bs = MakeSource[ Nx, Ny, dx, dy,RHS] ;
rhsvec = bx + by + bs;
D2mat = DX + DY;


Solving for the potential:


fxy = LinearSolve[D2mat, rhsvec];


This function gives any entry nm in a matrix an index that corresponds to the entry's position in a list.


g[n_,m_] := (m - 1)*Nx + n


The following constructs a table that asigns to every x and y point its corresponding potential.


Datatable= Table[{ii dx, jj dy, fxy[[ g[ii,jj] ]]},{ii,1, Nx}, {jj,1,Ny}];
Datatable = Partition[Flatten[Datatable],3];


Finally, with a table listing the values of the potential for every point (x,y), we plot the potential. 


Ga = ListContourPlot[Datatable,Contours-> 20,ColorFunction-> "GrayTones"]	


The plot of the potential is shown due to a point charge. The lines indicate the equipotentials, that is, where the potential is constant.


Part b


Nx=49;
Ny=50;
dx = 1/(Nx + 1);
dy = 1/(Ny+1);
RHS[x_,y_] :=Module[{ret},
ret = 0.0;
If[Round[x/dx]==25 && Round[y/dy]==25,ret=1];
If[Round[x/dx]==25 && Round[y/dy]==26,ret=-1];
Return[ret];
]


{DX, bx} =DXmat[Nx, Ny,dx,dy,  left, right] ;
{DY, by} = DYmat[Nx, Ny,  dx,dy,top, bottom];
bs = MakeSource[ Nx, Ny, dx, dy,RHS] ;
rhsvec = bx + by + bs;
D2mat = DX + DY;


fxy = LinearSolve[D2mat, rhsvec];


g[n_,m_] := (m - 1)*Nx + n


Datatable= Table[{ii dx, jj dy, fxy[[ g[ii,jj] ]]},{ii,1, Nx}, {jj,1,Ny}];
Datatable = Partition[Flatten[Datatable],3];


Ga = ListContourPlot[Datatable,Contours-> 20,ColorFunction-> "BlueGreenYellow"]	


Part c


Nx=50;
Ny=50;
dx = 1/(Nx + 1);
dy = 1/(Ny+1);
RHS[x_,y_] :=Module[{ret},
ret = 0.0;
If[Round[x/dx]==25 && Round[y/dy]==25,ret=1];
If[Round[x/dx]==26 && Round[y/dy]==26,ret=1];
If[Round[x/dx]==25 && Round[y/dy]==26,ret=-1];
If[Round[x/dx]==26 && Round[y/dy]==25,ret=-1];
Return[ret];
];
{DX, bx} =DXmat[Nx, Ny,dx,dy,  left, right] ;
{DY, by} = DYmat[Nx, Ny,  dx,dy,top, bottom];
bs = MakeSource[ Nx, Ny, dx, dy,RHS] ;
rhsvec = bx + by + bs;
D2mat = DX + DY;
fxy = LinearSolve[D2mat, rhsvec];


g[n_,m_] := (m - 1)*Nx + n
Datatable= Table[{ii dx, jj dy, fxy[[ g[ii,jj] ]]},{ii,1, Nx}, {jj,1,Ny}];
Datatable = Partition[Flatten[Datatable],3];
Ga = ListContourPlot[Datatable,Contours-> 20,ColorFunction-> "Rainbow"]	


Problem 4.12


This function will generate a box that will have charge.


Nx=60;
Ny=60;
dx = 1/(Nx + 1);
dy = 1/(Ny+1);
RHS[x_,y_] :=Module[{ret,i,j},
ret = 0.0;
i=15;
j=15;
While[i<=46,
If[Round[x/dx]==i && Round[y/dy]==j,ret=1];
i=i+1];
i=15;
j=46;
While[i<=46,
If[Round[x/dx]==i && Round[y/dy]==j,ret=1];
i=i+1];
i=15;
j=15;
While[j<=46,
If[Round[x/dx]==i && Round[y/dy]==j,ret=1];
j=j+1];
i=46;
j=15;
While[j<=46,
If[Round[x/dx]==i && Round[y/dy]==j,ret=1];
j=j+1];
Return[ret];
];

{DX, bx} =DXmat[Nx, Ny,dx,dy,  left, right] ;
{DY, by} = DYmat[Nx, Ny,  dx,dy,top, bottom];
bs = MakeSource[ Nx, Ny, dx, dy,RHS] ;
rhsvec = bx + by + bs;
D2mat = DX + DY;
fxy = LinearSolve[D2mat, rhsvec];


g[n_,m_] := (m - 1)*Nx + n
Datatable= Table[{ii dx, jj dy, fxy[[ g[ii,jj] ]]},{ii,1, Nx}, {jj,1,Ny}];
Datatable = Partition[Flatten[Datatable],3];
Ga = ListContourPlot[Datatable,Contours-> 20,ColorFunction-> "GrayTones"]	


Problem 4.14


This function will generate three parallel finite lines of charge. The top and the bottom ones will be positively charged, while the one in the middle will be negatively charged.


Nx=49;
Ny=49;
dx = 1/(Nx + 1);
dy = 1/(Ny+1);
RHS[x_,y_] :=Module[{ret,i,j},
ret = 0.0;
i=5;
j=25;
While[i<=46,
If[Round[x/dx]==i && Round[y/dy]==j,ret=-2];
i=i+1
];
i=5;
j=15;
While[i<=46,
If[Round[x/dx]==i && Round[y/dy]==j,ret=1];
i=i+1
];
i=5;
j=35;
While[i<=46,
If[Round[x/dx]==i && Round[y/dy]==j,ret=1];
i=i+1
];
Return[ret];
];

{DX, bx} =DXmat[Nx, Ny,dx,dy,  left, right] ;
{DY, by} = DYmat[Nx, Ny,  dx,dy,top, bottom];
bs = MakeSource[ Nx, Ny, dx, dy,RHS] ;
rhsvec = bx + by + bs;
D2mat = DX + DY;
fxy = LinearSolve[D2mat, rhsvec];


g[n_,m_] := (m - 1)*Nx + n
Datatable= Table[{ii dx, jj dy, fxy[[ g[ii,jj] ]]},{ii,1, Nx}, {jj,1,Ny}];
Datatable = Partition[Flatten[Datatable],3];
Ga = ListContourPlot[Datatable,Contours-> 20,ColorFunction-> "GrayTones"]	


This function will generate the gradient of a potential by approximating the values of the derivative of the potential. 


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


We can then obtain the electric field.


Efield = N[GradV[fxy,Nx, Ny,dx,dy,-1]];
GE = ListVectorPlot[Efield, VectorStyle->GrayLevel[0.5]]


As expected, the electric field in-between the top an bottom plates is non-zero, while the electric field above and below the top and bottom plates is zero. We can also see that the electric fields curves at the edges of the plates.


Problem 4.16


For this problem, in the functions the variable x represents s, and the variable y represents \[Phi].


The next three functions were adapted to work for polar coordinates. For the MakeSource function the function g[n,m] was changed to an appropiate one for polar coordiantes. For both of the derivative functions extra terms in the coefficients of the matrix were added, following the formula for the gradient for polar coordinates.


MakeSource[ Nx_, Ny_, dx_, dy_, source_] := Module[{svec, g,j, index, indey},
g[n_,m_] := Mod[m-1,Ny]Nx+n;

svec = Table[0.0,{i,1,Nx Ny}];
For[index = 1, index <= Nx, index=index+1,
For[indey = 1, indey <= Ny, indey = indey + 1,
j = g[index,indey];
svec[[j ]] = source[ index dx, indey dy];
];
];
Return[svec];
]


DSmat[Nx_, Ny_, dx_, dy_, left_, right_] := Module[{g,retmat,bside, index, indey,i,j,k},
g[n_,m_] := Mod[m-1,Ny]Nx+n;

retmat = Table[0.0,{i,1,Nx Ny},{j,1,Nx Ny}];
bside = Table[0.0,{i,1,Nx Ny}];

For[indey = 1, indey <= Ny, indey = indey+1,
j = g[1,indey];
k = g[2,indey];

retmat[[j,j ]] += -2.0/dx^2;
retmat[[j,k ]] += 1.0/dx^2 +1/(2 dx( s0+1 dx));
bside[[ j ]] += -left[indey dy](1/dx^2 -1/(2 dx (s0+1 dx)));

For[index = 2, index <= Nx-1, index = index+1,
i = g[index-1,indey];
j = g[index,indey];
k = g[index + 1,indey];
retmat[[j,i]] += 1.0/dx^2-1/(2 dx (s0+index dx));
retmat[[j,j]] += -2.0/dx^2;
retmat[[j,k]] += 1.0/dx^2+1/(2 dx (s0+index dx));
];
i = g[Nx-1,indey];
j = g[Nx,indey];
retmat[[j,i ]] += 1.0/dx^2-1/(2 dx (s0+Nx dx));
retmat[[ j,j]] += -2.0/dx^2;
bside[[ j]] += -right[indey dy](1/dx^2+1/(2 dx (s0+Nx dx)));
];
Return[{retmat, bside}];
]


Dphimat[Nx_, Ny_, dx_, dy_] := Module[{g,retmat,bside, index, indey,i,j,k},
g[n_,m_] := Mod[m-1,Ny]Nx+n;

retmat = Table[0.0,{i,1,Nx Ny},{j,1,Nx Ny}];
bside = Table[0.0,{i,1,Nx Ny}];

For[index = 1, index <= Nx, index = index+1,
j = g[index,1];
k = g[index, 2];

For[indey = 1, indey <= Ny, indey = indey+1,
i = g[index,indey-1];
j = g[index,indey];
k = g[index ,indey+1];
retmat[[j,i]] += 1.0/(dy^2 ( (s0+index dx))^2);
retmat[[j,j]] += -2.0/(dy^2((s0+index dx))^2);
retmat[[j,k]] += 1.0/(dy^2((s0+index dx))^2);
];

];
Return[{retmat, bside}];
]



We will obtain the potential outside a sphere. The sphere is at a constant potential Subscript[V, o]=10.


RHS[x_,y_] := 0;
top[x_] :=0.;
bottom[x_] :=0.;
left[y_] := 10.;
right[y_] :=0.


The sphere will have radius 0.1. 


s0=0.1;
sf=2;
Ny=40;
Nx=20;
dx = (sf-s0)/(Nx+1);
dy =  2Pi/Ny;


{DX, bx} = DSmat[Nx, Ny, dx, dy, left, right] ;
{DY, by} = Dphimat[Nx, Ny, dx, dy];
bs = MakeSource[ Nx, Ny, dx, dy,RHS] ;
rhsvec = bx + by + bs;
D2mat = DX + DY;
fxy = LinearSolve[D2mat, rhsvec];


The following functions were also adapted to work with polar coordinates. In the table, the first entry is the radius, which is the cosine of the initial radius s0 plus the distance s. The second entry is the sine of the initial radius s0 plus the distance s.


g[n_,m_] := Mod[m-1,Ny]Nx+n;
Datatable= Table[{(s0 +ii dx) Cos[jj dy],(s0 +ii dx) Sin[jj dy], fxy[[ g[ii,jj] ]]},{ii,1, Nx}, {jj,1,Ny}];
Datatable = Partition[Flatten[Datatable],3];

Ga = ListContourPlot[Datatable,Contours-> 20,ColorFunction-> "GrayTones"]


This table will define the potential in function of s for \[Phi] = 0:


Pot=Table[{(s0 +ii dx) , fxy[[ g[ii,1] ]]},{ii,1, Nx}];
Datatable = Partition[Flatten[Datatable],2];
DPot=ListPlot[Pot];


This is the function obtained theoretically:


TPot[s_]:=10 Log[s/2]/Log[0.1/2];
ThPot=Plot[TPot[s],{s,0,2}];


Both the potential as a function of s obtained by numerical methods and analytically are shown in the plot below. 


Show[DPot,ThPot]


Problem 4.17


This function will generate a non-uniform potential on the surface of our sphere. This potential will be 10 Sin [\[Phi]]. 


RHS[x_,y_] := 0;
top[x_] :=0.;
bottom[x_] :=0.;
left[y_] := 10 Sin[y];
right[y_] :=0.


s0=0.1;
sf=2;
Ny=40;
Nx=20;
dx = (sf-s0)/(Nx+1);
dy =  2Pi/Ny;


{DX, bx} = DSmat[Nx, Ny, dx, dy, left, right] ;
{DY, by} = Dphimat[Nx, Ny, dx, dy];
bs = MakeSource[ Nx, Ny, dx, dy,RHS] ;
rhsvec = bx + by + bs;
D2mat = DX + DY;
fxy = LinearSolve[D2mat, rhsvec];


g[n_,m_] := Mod[m-1,Ny]Nx+n;
Datatable= Table[{(s0 +ii dx) Cos[jj dy],(s0 +ii dx) Sin[jj dy], fxy[[ g[ii,jj] ]]},{ii,1, Nx}, {jj,1,Ny}];
Datatable = Partition[Flatten[Datatable],3];

Ga = ListContourPlot[Datatable,Contours-> 20,ColorFunction-> "GrayTones"]
