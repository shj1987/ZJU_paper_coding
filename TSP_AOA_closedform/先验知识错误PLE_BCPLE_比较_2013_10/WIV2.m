function [P_Bias, H_Bias] = WIV2(SB,DOA,UL,NA,UU)  %%%%%伪线性最小二乘法WIV
% DOA : 1 x n, t时刻各个阵列接收到绝对坐标系下的角度信息;
% Array : n x 2， t时刻各个阵列的位置;
%UBC是bias compensated PLE
lt=length(SB(:,1));   
for j=1:lt
    BCN=SB(j,:); %beacon
    aj=BCN(1);
    bj=BCN(2);
    Qj=[aj+bj*UU(3),bj-aj*UU(3)];
    AOA(j)= atan2(UU(2)-Qj(2), UU(1)-Qj(1));
    W(j,j)=(Qj(1)-UU(1))^2+(Qj(2)-UU(2))^2;
end
for j=1:lt
    L(j,1) = sin(AOA(j));
    L(j,2) = -cos(AOA(j));
    L(j,3) = -(SB(j,1)*cos(AOA(j))+SB(j,2)*sin(AOA(j)));
end
for j=1:lt
    C(j,1) = sin(DOA(j));
    C(j,2) = -cos(DOA(j));
    C(j,3) = -(SB(j,1)*cos(DOA(j))+SB(j,2)*sin(DOA(j)));
    D(j,1) = sin(DOA(j))*SB(j,1) - cos(DOA(j))*SB(j,2);
end
IW=inv(W);
T = inv(L'*IW*C)*L'*IW*D;
X=(T(1)-T(3)*T(2))/(1+T(3)^2);
Y=(T(2)+T(1)*T(3))/(1+T(3)^2);
Fai=atan(T(3));
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

