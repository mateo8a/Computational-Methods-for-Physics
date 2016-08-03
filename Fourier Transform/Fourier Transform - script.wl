(* ::Package:: *)

Lab 7


Mateo Ochoa


Problem 7.12


DiscreteFourier[list_]:=Module[{Plist, k,Ns},
Ns=Length[list];
Plist = Table[0.,{i,1,Ns}];
For[k=0,k<=Ns-1,k=k+1,
Plist[[k+1]] =Sum[Exp[2 Pi I j k/Ns] list[[j+1]],{j,0,Ns-1}]; 
];
Return[Plist];
]


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


p[t_]:=Sin[2 Pi 5 t];
n=2^8;
dt=2/n;
df=0.5;


Using the discrete Fourier transform:


data=Table[p[t],{t,0.,dt(n-1),dt}];
ftransform=DiscreteFourier[data];
ftransform=RotateRight[ftransform,n/2];
pspec=Table[{k df,(Abs[ftransform[[(k+n/2)+1]]])^2},{k,-n/2,n/2 -1}];


ListLinePlot[pspec,PlotRange->{{-10,10},All},AxesLabel->{"\!\(\*SubscriptBox[\(f\), \(k\)]\)","|P(\!\(\*SubscriptBox[\(f\), \(k\)]\))\!\(\*SuperscriptBox[\(|\), \(2\)]\)"}]


Using the fast Fourier transform:


FTTtransform=FFT[data];
FTTtransform=RotateRight[FTTtransform,n/2];
pspec=Table[{k df,(Abs[FTTtransform[[(k+n/2)+1]]])^2},{k,-n/2,n/2 -1}];


ListLinePlot[pspec,PlotRange->{{-10,10},All},AxesLabel->{"\!\(\*SubscriptBox[\(f\), \(k\)]\)","|P(\!\(\*SubscriptBox[\(f\), \(k\)]\))\!\(\*SuperscriptBox[\(|\), \(2\)]\)"}]


Replacing the second half of the data with zeroes, we have:


data=Table[If[t<=dt (n-1)/2,p[t],0],{t,0.,dt(n-1),dt}];
ftransform=DiscreteFourier[data];
ftransform=RotateRight[ftransform,n/2];
pspec=Table[{k df,(Abs[ftransform[[(k+n/2)+1]]])^2},{k,-n/2,n/2 -1}];


ListLinePlot[pspec,PlotRange->{{-10,10},All},AxesLabel->{"\!\(\*SubscriptBox[\(f\), \(k\)]\)","|P(\!\(\*SubscriptBox[\(f\), \(k\)]\))\!\(\*SuperscriptBox[\(|\), \(2\)]\)"}]


Problem 7.13


data=Import["http://people.reed.edu/~jfrankli/Book/downloads/files/ch7.dat"];


psubj=Table[data[[i,2]],{i,1,Length[data]}];


n=Length[psubj];
dt=data[[2,1]]-data[[1,1]];
df = 1/(Length[psubj] dt);
FTTtransform=FFT[psubj];
FTTtransform=RotateRight[FTTtransform,n/2];
pspec=Table[{k df,(Abs[FTTtransform[[(k+n/2)+1]]])^2},{k,-n/2,n/2 -1}];


ListLinePlot[pspec,PlotRange->{{-30,30},{0,100000}},AxesLabel->{"\!\(\*SubscriptBox[\(f\), \(k\)]\)","|P(\!\(\*SubscriptBox[\(f\), \(k\)]\))\!\(\*SuperscriptBox[\(|\), \(2\)]\)"}]


We obtain the locations (in the grid) of the relevant frequencies, and the magnitudes of the Fourier transform at those frequencies with:


Fk[data_]:=Module[{i,list,listf,k,kk},
list=Table[0,{i,1,Length[data]}];
k=1;
For[i=1,i<=Length[data],i=i+1,
If[(Abs[data[[i]]])^2>1000,list[[k]]={i-1,data[[i]]};k=k+1]
];
kk=k-1;
listf=Table[0,{i,1,kk}];
For[i=1,i<=kk,i=i+1,
listf[[i]]=list[[i]]];
Return[listf]]


fk=Fk[FTTtransform];


To obtain the frequencies in Hz, we implement the function:


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


The frequencies and the coefficient (Fourier transform) corresponding to each frequency of the signal will be:


fnc=Chop[freqHz[fk]];
MatrixForm[fnc]


We can know which frequencies correspond to sine or cosine functions by looking if it is real or imaginary.


The coefficients for the Cosine and Sine functions will be given by:


Coeff[data_]:=Module[{transform,list,A1},
list=Table[0,{i,1,Length[data]}];
transform=Table[data[[i,2]],{i,1,Length[data]}];
For[i=1,i<=Length[transform],i=i+1,
If[Element[transform[[i]],Reals] ,A1=2*transform[[i]]*dt*df,
A1=-2*I*transform[[i]]*dt*df;
];
list[[i]]={data[[i,1]],A1}];
Return[Chop[list]];
]


coefA=Coeff[fnc];
MatrixForm[coefA]


Here we have every frequency with the coefficient that will go in front of the sine or cosine function.


The original function can be constructed by considering only the positive frequencies.


originalf[t_]:=coefA[[8,2]] Cos[2 Pi coefA[[8,1]] t]+
coefA[[5,2]] Sin[2 Pi coefA[[5,1]] t]+
coefA[[6,2]] Sin[2 Pi coefA[[6,1]] t]+
coefA[[7,2]] Sin[2 Pi coefA[[7,1]] t]


Plotting the constructed function from 0 to 2:


Plot[originalf[t],{t,0,2},AxesLabel->{"t","f(t)"}]


Compared to the original data.


ListLinePlot[data,AxesLabel->{"t","f(t)"}]


Problem 7.15


The inverse FFT will be:


inverseFFT[invec_,NN_]:=Module[{iFFT,Den},
iFFT[invec] := Module[{k,elist,olist,netlist,eft,oft,n,omega},
If[Length[invec] ==  1, 
Return[invec];
];
n = Length[invec];
omega = 1;
elist = Table[invec[[j]],{j,1,n,2}];
olist = Table[invec[[j]],{j,2,n,2}];
eft = inverseFFT[elist,NN];
oft = inverseFFT[olist,NN];
netlist = Table[0,{j,0,n-1}];
For[k = 0, k <= n/2 - 1, k = k+1,
netlist[[k+1]] = eft[[k+1]] + omega oft[[k+1]];
netlist[[k + n/2+1]] = eft[[k+1]] - omega oft[[k+1]];
omega = omega Exp[-2 Pi I/n];
];
Return[netlist];
];
Den=If[Length[invec]<NN,1,NN];
Return[iFFT[invec]/Den];
]


I modified the original FTT by adding the negative sign for the omega expression. For the division, I define an input number which will always be the length of invec. Since the function iFFT calls itself inside inverseFFT, the length of invec changes every iteration. NN will save the original length regarless of how many iterations iFFT goes through.


We will download data from the course site and FTT it.


data=Import["http://people.reed.edu/~jfrankli/Book/downloads/files/ch7.dat"];


psubj=Table[data[[i,2]],{i,1,Length[data]}];
FTTtransform=FFT[psubj];
dt=data[[2,1]]-data[[1,1]];


Our inverseFFT will give:


inverse=Chop[inverseFFT[FTTtransform,Length[FTTtransform]]];
ListLinePlot[Table[{i dt,inverse[[i]]},{i,1,Length[inverse]}],AxesLabel->{"t","f(t)"}]


Plotting the original data we have:


ListLinePlot[data,AxesLabel->{"t","f(t)"}]


Problem 7.16


The band pass filter will be given by:


BandPass[invec_,cutbelow_,cutabove_, T_] := Module[{outvec,NN,dt,df,cutbindex,cutaindex,FTdata,index},
NN = Length[invec];
dt = N[T/NN];
df= N[1/(NN dt)];
FTdata = FFT[N[invec]];
cutbindex = Round[cutbelow/df];
cutaindex = Round[cutabove/df];
outvec = Table[0.0,{j,1,Length[FTdata]}];
For[index = cutbindex,index <= cutaindex, index=index+1,
outvec[[index]] = FTdata[[index]];
outvec[[Length[FTdata]-(index-2)]] = Conjugate[FTdata[[index]]];
];
outvec = inverseFFT[outvec,Length[outvec]];
Return[Chop[outvec]];
]


The index for which the For loop will start will be given by the value of cutbelow. 


data = Import["http://people.reed.edu/~jfrankli/Book/downloads/files/GraviteaTime.wav","Data"][[1]];


This function will find the most appropiate n for the exponent in 2^n.


nfinder[data_]:=Module[{n},
n=1;
While[2^n<Length[data],
n=n+1];
Return[n-1];]


sample=Table[data[[j]],{j,1,2^nfinder[data]}];
dt = 1/41000;
df = 1/(Length[sample] dt);


We used the BandPass on the imported data.


filtereddata=BandPass[sample,440,880,Length[sample] dt];


ListPlay[Chop[filtereddata],SampleRate->41000]


n=Length[filtereddata];
FTTtransform=FFT[filtereddata];
FTTtransform=RotateRight[FTTtransform,n/2];
pspec=Table[{k df,(Abs[FTTtransform[[(k+n/2)+1]]])^2},{k,-n/2,n/2 -1}];


ListLinePlot[pspec,PlotRange->{{0,2000},{0,10^7}},AxesLabel->{"\!\(\*SubscriptBox[\(f\), \(k\)]\)","|P(\!\(\*SubscriptBox[\(f\), \(k\)]\))\!\(\*SuperscriptBox[\(|\), \(2\)]\)"}]
