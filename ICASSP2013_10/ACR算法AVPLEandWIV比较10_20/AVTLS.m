function [OV2,P_Bias, H_Bias] = AVTLS(SB,DOA,UL,NA)  %%%%%Total Least Square伪线性最小二乘法

% DOA : 1 x n, t时刻各个阵列接收到绝对坐标系下的角度信息;
% Array : n x 2， t时刻各个阵列的位置;
lt=length(DOA);          
for i=1:lt
    C(i,1) = sin(DOA(i));
    C(i,2) = -cos(DOA(i));
    C(i,3) = -(SB(i,1)*cos(DOA(i))+SB(i,2)*sin(DOA(i)));
    D(i,1) = sin(DOA(i))*SB(i,1) - cos(DOA(i))*SB(i,2);
end
CD=[C,D];       %增广矩阵
[U,S,V]=svd(CD);
V4=V(:,4);
T = -V4(1:3)/V4(4);
% TC=CD'*CD;
% Val=eig(TC'*TC);
% sgm=min(Val);    %最小特征值
% T2 = inv(C'*C-sgm*eye(3))*C'*D
X=(T(1)-T(3)*T(2))/(1+T(3)^2);
Y=(T(2)+T(1)*T(3))/(1+T(3)^2);
Fai=atan(T(3));
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
P_Bias=sqrt((UL(1)-X)^2+(UL(2)-Y)^2);
H_Bias=abs(Fi-NA)*180/pi;
OV2=[X,Y,Fi];
