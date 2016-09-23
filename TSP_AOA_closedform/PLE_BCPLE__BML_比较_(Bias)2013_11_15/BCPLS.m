function [UBC,OV,PH] = BCPLS(SB,DOA,UL,NA)  %%%%%伪线性最小二乘法PLS
global SIGMA
sgm=(SIGMA)*pi/180;
lt=length(DOA);          
for i=1:lt
    C(i,1) = sin(DOA(i));
    C(i,2) = -cos(DOA(i));
    C(i,3) = -(SB(i,1)*cos(DOA(i))+SB(i,2)*sin(DOA(i)));
    D(i,1) = sin(DOA(i))*SB(i,1) - cos(DOA(i))*SB(i,2);
end
T = inv(C'*C)*C'*D;
AT1=0;
AT2=0;
AT3=0;
%%%%%aj,bj denote the location of beacon%%%%%%%%
for j=1:lt
    BN=SB(j,:);
    aj=BN(1);
    bj=BN(2);
    cj=aj^2+bj^2;
    AT1=(T(1)-T(3)*bj-aj)+AT1;
    AT2=(T(2)-bj+T(3)*aj)+AT2;
    AT3=T(2)*aj-T(1)*bj+T(3)*cj+AT3;
end
C11=0;
C12=0;
C13=0;
C21=0;
C22=0;
C23=0;
C31=0;
C32=0;
C33=0;
for ni=1:lt
    SBB=SB(ni,:);
    C11=C11+sin(DOA(ni))^2+sgm^2*cos(2*DOA(ni));
    C12=C12-0.5*sin(2*DOA(ni))+sgm^2*sin(2*DOA(ni));
    C13=C13-0.5*SBB(1)*sin(2*DOA(ni))-SBB(2)*sin(DOA(ni))^2+sgm^2*(SBB(1)*sin(2*DOA(ni))-SBB(2)*cos(2*DOA(ni)));
    C21=C21-0.5*sin(2*DOA(ni))+sgm^2*sin(2*DOA(ni));
    C22=C22+cos(DOA(ni))^2-sgm^2*cos(2*DOA(ni));
    C23=C23+SBB(1)*cos(DOA(ni))^2+0.5*SBB(2)*sin(2*DOA(ni))-sgm^2*(SBB(1)*cos(2*DOA(ni))+SBB(2)*sin(2*DOA(ni)));
    C31=C31-0.5*SBB(1)*sin(2*DOA(ni))-SBB(2)*sin(DOA(ni))^2+sgm^2*(SBB(1)*sin(2*DOA(ni))-SBB(2)*cos(2*DOA(ni)));
    C32=C32+SBB(1)*cos(DOA(ni))^2+0.5*SBB(2)*sin(2*DOA(ni))-sgm^2*(SBB(1)*cos(2*DOA(ni))+SBB(2)*sin(2*DOA(ni)));
    C33=C33+(SBB(1)*cos(DOA(ni))+SBB(2)*sin(DOA(ni)))^2+sgm^2*(cos(2*DOA(ni))*(SBB(2)^2-SBB(1)^2)-2*SBB(1)*SBB(2)*sin(2*DOA(ni)));
end
EC=[C11 C12 C13;C21 C22 C23;C31 C32 C33];
AT=[AT1,AT2,AT3]';
%UBC=T+inv(C'*C)*AT*sgm^2;
UBC=T+inv(EC)*AT*sgm^2;  
X=(UBC(1)-UBC(3)*UBC(2))/(1+UBC(3)^2);
Y=(UBC(2)+UBC(1)*UBC(3))/(1+UBC(3)^2);
Fai=atan(UBC(3)); 
Tha=atan2(SB(1,2)-Y,SB(1,1)-X);
if  Fai<0
    Fi=Fai+pi;
    if rem(abs(Tha-DOA(1)-Fi),2*pi)>pi/2
        Fi=Fai;
    end
else 
     Fi=Fai;
     if rem(abs(Tha-DOA(1)-Fi),2*pi)>pi/2
        Fi=Fai+pi;
    end
end
P_Bias=sqrt((UL(1)-X)^2+(UL(2)-Y)^2);
H_Bias=abs(Fi-NA)*180/pi;
OV=[X,Y,Fai];
PH=[X,Y,Fai];
