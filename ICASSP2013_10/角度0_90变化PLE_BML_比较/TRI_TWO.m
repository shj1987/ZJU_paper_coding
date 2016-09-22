function [P_Bias, H_Bias] = TRI_TWO(SB,DOA,UL,NA)  %%%%%伪线性最小二乘法PLS
% DOA : 1 x n, t时刻各个阵列接收到绝对坐标系下的角度信息;
% Array : n x 2， t时刻各个阵列的位置;
% global SIGMA
lt=length(DOA);  
VT=nchoosek(1:lt,2);
ko=0;
L=length(VT(:,1));
for j=1:L
    VC=VT(j,:);
    j1=VC(1);
    j2=VC(2);
    BN1=SB(j1,:);
    BN2=SB(j2,:);
    AOA1=abs(DOA(j1)-DOA(j2));
    if AOA1>pi
        AOA1=2*pi-AOA1;
    end
    if AOA1> pi/2
        AOA1=2*pi-2*AOA1;
    end
    AOA=2*AOA1;
    d2=(BN1(1)-BN2(1))^2+(BN1(2)-BN2(2))^2;
    r2=d2/(2-2*cos(AOA));
    q=0.5*(BN1(1)^2-BN2(1)^2+BN1(2)^2-BN2(2)^2)/(BN2(1)-BN1(1));
    p=(BN2(2)-BN1(2))/(BN2(1)-BN1(1));
    ww=2*(q*p+p*BN1(1)-BN1(2));
    z=(q+BN1(1))^2+BN1(2)^2-r2;
    y1=0.5*(-ww+sqrt(ww^2-4*z*(p^2+1)))/(p^2+1);
    y2=0.5*(-ww-sqrt(ww^2-4*z*(p^2+1)))/(p^2+1);
    x1=-p*y1-q;
    x2=-p*y2-q;
    O1=[];
    O1=[x1,y1;x1,y2;x2,y1;x2,y2];
%     O1=unique(xy,'rows');%排序
    dm=150;
    ko=ko+1;
    for ii=1:length(O1(:,1))
        OO=O1(ii,:);
        Xr=sqrt((OO(1)-UL(1))^2+(OO(2)-UL(2))^2);
        dd=abs(Xr-sqrt(r2));%radius
        if dd<dm
           dm=dd;
           O(ko,:)=OO;
           DR(ko)=r2;  
        end
    end
end
XE=[];
RN=[];
% XN=O; % last one beacon
% RN=DR;  %r^2
kl=length(DR);
XN=O(kl,:); % last one beacon
RN=DR(kl);  %r^2
for ij=1:kl-1
    XB=O(ij,:); %beacon
    C(ij,1) = 2*(XB(1)-XN(1));
    C(ij,2) = 2*(XB(2)-XN(2));
    D(ij,1) = XB(1)^2-XN(1)^2+XB(2)^2-XN(2)^2-DR(ij)+RN;
end
XY = inv(C'*C)*C'*D;
ADOA=atan2(SB(:,2)-XY(2),SB(:,1)-XY(1));
kd=0;
for di=1:length(ADOA)
    ORI=ADOA(di);
    Hd=ORI-DOA(di);
    kd=kd+1;
    HB(kd)=abs(Hd-NA)*180/pi
end
%H_Bias=abs(Hd-NA)*180/pi;
H_Bias=mean(HB);
P_Bias=sqrt((UL(1)-XY(1))^2+(UL(2)-XY(2))^2);