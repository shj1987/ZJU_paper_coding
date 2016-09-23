function [OV,P_Bias, H_Bias] = PLS(SB,DOA,UL,NA)  %%%%%伪线性最小二乘法PLS

% DOA : 1 x n, t时刻各个阵列接收到绝对坐标系下的角度信息;
% Array : n x 2， t时刻各个阵列的位置;
lt=length(DOA);          
for i=1:lt
    C(i,1) = sin(DOA(i));
    C(i,2) = -cos(DOA(i));
    C(i,3) = -(SB(i,1)*cos(DOA(i))+SB(i,2)*sin(DOA(i)));
    D(i,1) = sin(DOA(i))*SB(i,1) - cos(DOA(i))*SB(i,2);
end
T = inv(C'*C)*C'*D;
X=(T(1)-T(3)*T(2))/(1+T(3)^2);
Y=(T(2)+T(1)*T(3))/(1+T(3)^2);
Fai=atan(T(3)); 
if Fai<0
    Fai=Fai+pi;            %修正角度
end
P_Bias=sqrt((UL(1)-X)^2+(UL(2)-Y)^2);
H_Bias=abs(Fai-NA)*180/pi;
OV=[X,Y,Fai];