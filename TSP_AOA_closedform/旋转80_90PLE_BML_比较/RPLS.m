function [T,P_B, H_B] = RPLS(SB,DOA,UL,NA)  %%%%%伪线性最小二乘法PLS
% DOA : 1 x n, t时刻各个阵列接收到绝对坐标系下的角度信息;
% Array : n x 2， t时刻各个阵列的位置;
lt=length(DOA);          
for i=1:lt
    C(i,1) = sin(DOA(i));
    C(i,2) = -cos(DOA(i));
    C(i,3) = -SB(i,1)*sin(DOA(i))+SB(i,2)*cos(DOA(i));
    D(i,1) = -cos(DOA(i))*SB(i,1) - sin(DOA(i))*SB(i,2);
end
T = inv(C'*C)*C'*D;
X=(T(2)+T(3)*T(1))/(1+T(3)^2);
Y=(T(2)*T(3)-T(1))/(1+T(3)^2);
Fai=atan(T(3))+pi/2;
% Tha=atan2(SB(1,2)-Y,SB(1,1)-X);
% if Tha<0
%    Tha= Tha+2*pi;
% end
% if  Fai<0
%     Fi=Fai+pi;
%     if abs(Tha-DOA(1)-Fi)>pi/2
%         Fi=Fai;
%     end
% else 
%      Fi=Fai;
%       if abs(Tha-DOA(1)-Fi)>pi/2
%         Fi=Fai+pi;
%     end
% end
P_B=sqrt((UL(1)-X)^2+(UL(2)-Y)^2);
H_B=abs(Fai-NA)*180/pi;
