function [flg,T,P_Bias, H_Bias] = PLS(SB,DOA,UL,NA)  %%%%%伪线性最小二乘法PLS

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
for i=1:lt
   Tha=atan2(SB(i,2)-Y,SB(i,1)-X);
   if Tha<0
      Tha= Tha+2*pi;
   end
   Ag(i)=Tha;
end
zg=0;
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
cnt=0;
thred=pi/30; %与cnt 相反
for i=1:lt
    da=abs(Ag(i)-DOA(i)-Fi);
    if da>thred
       cnt=cnt+1;
    end
end

if cnt>=3
    flg=1;
else
    flg=0;
end
%%40-7/ 35-4
P_Bias=sqrt((UL(1)-X)^2+(UL(2)-Y)^2);
H_Bias=abs(Fi-NA)*180/pi;

