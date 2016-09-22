%true emitter location
p=[55;35]

%receiver path parameters y = a*x + b
a=-.2;
b=14;

%no. of bearing measurements
N=20;

%observer x-axis range
rxl=5;
rxu=45;

%geometry plot axis ranges
xy=[0, 60, 0, 40];

%CTLS parameters
mrec=100; %maximum iterations for successive projections
errthr=1e-4; %error threshold for convergence


%ML parameters
pinit=[45;25];

%Noise and parameter data
%bearing std, observer pos std, stepsize, no of iterations
         
par=[1*pi/180,.1,9e-4,100;
    2*pi/180,.1,3e-3,100;
    4*pi/180,.1,1e-2,100;
    6*pi/180,.1,2e-2,100;
    8*pi/180,.1,4e-2,100;
    2*pi/180,0,4e-3,100;
    2*pi/180,.2,4e-3,100;
    2*pi/180,.4,4e-3,100;
    2*pi/180,.6,4e-3,100;
    2*pi/180,.8,4e-3,100;
    2*pi/180,1.0,4e-3,100];
