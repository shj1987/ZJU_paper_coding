function [P_B, H_B] = RAVTLS(SB,DOA,UL,NA)  %%%%%伪线性最小二乘法PLS
% DOA : 1 x n, t时刻各个阵列接收到绝对坐标系下的角度信息;
% Array : n x 2， t时刻各个阵列的位置;
lt=length(DOA);     
for i=1:lt
    C(i,1) = sin(DOA(i));
    C(i,2) = -cos(DOA(i));
    C(i,3) = -SB(i,1)*sin(DOA(i))+SB(i,2)*cos(DOA(i));
    D(i,1) = -cos(DOA(i))*SB(i,1) - sin(DOA(i))*SB(i,2);
end
CD=[C,D];       %增广矩阵
TC=CD'*CD;
Val=eig(TC'*TC);
sgm=min(Val);    %最小特征值
T = inv(C'*C-sgm*eye(3))*C'*D; 
X=(T(2)+T(3)*T(1))/(1+T(3)^2);
Y=(T(2)*T(3)-T(1))/(1+T(3)^2);
Fai=atan(T(3))+pi/2;
for i=1:lt
   Tha=atan2(SB(i,2)-Y,SB(i,1)-X);
   if Tha<0
      Tha= Tha+2*pi;
   end
   Ag(i)=Tha;
end
if  Fai<0
    Fi=Fai+pi;
    if abs(Ag(1)-DOA(1)-Fi)>pi/2
        Fi=Fai; %0 degree
    end
else 
     Fi=Fai;
    if abs(Ag(1)-DOA(1)-Fi)>pi/2
        Fi=Fai+pi;
    end
end
P_B=sqrt((UL(1)-X)^2+(UL(2)-Y)^2);
H_B=abs(Fai-NA)*180/pi;
