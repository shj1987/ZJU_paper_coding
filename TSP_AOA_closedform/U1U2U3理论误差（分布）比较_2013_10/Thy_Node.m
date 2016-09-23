function EPLE=Thy_Node(TSB,DOA,Unode,TE)
global SIGMA
Sgm=SIGMA*pi/180;
W=length(DOA);
G1=0;
G2=0;
G3=0;
C11=0;
C12=0;
C13=0;
C21=0;
C22=0;
C23=0;
C31=0;
C32=0;
C33=0;
for ni=1:W
    SBB=TSB(ni,:);
    C11=C11+sin(DOA(ni))^2+Sgm^2*cos(2*DOA(ni));
    C12=C12-0.5*sin(2*DOA(ni))+Sgm^2*sin(2*DOA(ni));
    C13=C13-0.5*SBB(1)*sin(2*DOA(ni))-SBB(2)*sin(DOA(ni))^2+Sgm^2*(SBB(1)*sin(2*DOA(ni))-SBB(2)*cos(2*DOA(ni)));
    C21=C21-0.5*sin(2*DOA(ni))+Sgm^2*sin(2*DOA(ni));
    C22=C22+cos(DOA(ni))^2-Sgm^2*cos(2*DOA(ni));
    C23=C23+SBB(1)*cos(DOA(ni))^2+0.5*SBB(2)*sin(2*DOA(ni))-Sgm^2*(SBB(1)*cos(2*DOA(ni))+SBB(2)*sin(2*DOA(ni)));
    C31=C31-0.5*SBB(1)*sin(2*DOA(ni))-SBB(2)*sin(DOA(ni))^2+Sgm^2*(SBB(1)*sin(2*DOA(ni))-SBB(2)*cos(2*DOA(ni)));
    C32=C32+SBB(1)*cos(DOA(ni))^2+0.5*SBB(2)*sin(2*DOA(ni))-Sgm^2*(SBB(1)*cos(2*DOA(ni))+SBB(2)*sin(2*DOA(ni)));
    C33=C33+(SBB(1)*cos(DOA(ni))+SBB(2)*sin(DOA(ni)))^2+Sgm^2*(cos(2*DOA(ni))*(SBB(2)^2-SBB(1)^2)-2*SBB(1)*SBB(2)*sin(2*DOA(ni)));
end
U1=Unode(1)+Unode(2)*tan(Unode(3));
U2=Unode(2)-Unode(1)*tan(Unode(3));
U3=tan(Unode(3));
for hi=1:W
    HP=TSB(hi,:);
    aj=HP(1);
    bj=HP(2);
    G1=G1+U1-bj*U3-aj;
    G2=G2+U2-bj+aj*U3;
    G3=G3+(aj^2+bj^2)*U3+aj*U2-bj*U1;
end
EC=[C11 C12 C13;C21 C22 C23;C31 C32 C33];
G=[G1;G2;G3];
EPLE=-inv(EC)*G*Sgm^2;                  %求出误差
TU=[U1 U2 U3];
Rea_ERR=TE-TU;
TR=TU+EPLE' ;                              %估计计算出来的准确值
TX=(TR(1)-TR(3)*TR(2))/(1+TR(3)^2);
TY=(TR(2)+TR(1)*TR(3))/(1+TR(3)^2);
Fai=atan(TR(3)); 
Tha=atan2(TSB(1,2)-TY,TSB(1,1)-TX);
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
ErrL=sqrt((Unode(1)-TX)^2+(Unode(2)-TY)^2);
ErrH=abs(Fi-Unode(3))*180/pi;



