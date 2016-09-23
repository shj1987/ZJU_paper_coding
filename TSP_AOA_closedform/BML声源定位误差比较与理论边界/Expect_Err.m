function Err=Expect_Err(Ag,Bcon,Source)
global SIGMA
Sgma=SIGMA*pi/180;
N=length(Ag);
A11=0;
A12=0;
A21=0;
A22=0;
for ni=1:N
    A11=A11+sin(Ag(ni))^2+Sgma^2*cos(2*Ag(ni));
    A12=A12-0.5*sin(2*Ag(ni))+Sgma^2*sin(2*Ag(ni));
    A21=A21-0.5*sin(2*Ag(ni))+Sgma^2*sin(2*Ag(ni));
    A22=A22+cos(Ag(ni))^2-Sgma^2*cos(2*Ag(ni));
end
EA=[A11 A12;A21 A22]/N;
BS=mean(Bcon)-Source;
EAA=inv(EA)*BS'*Sgma^2;
Err=sqrt(EAA(1)^2+EAA(2)^2);