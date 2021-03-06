function [UBC,OV,P_Bias, H_Bias] = RBCPLS(SB,DOA,UL,NA)  %%%%%伪线性最小二乘法PLS
global SIGMA
sgm=(SIGMA)*pi/180;
lt=length(DOA);                
for i=1:lt
    C(i,1) = sin(DOA(i));
    C(i,2) = -cos(DOA(i));
    C(i,3) = -SB(i,1)*sin(DOA(i))+SB(i,2)*cos(DOA(i));
    D(i,1) = -cos(DOA(i))*SB(i,1) - sin(DOA(i))*SB(i,2);
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
    AT1=(T(1)-T(3)*aj+bj)+AT1;
    AT2=(T(2)-aj-T(3)*bj)+AT2;
    AT3=-T(1)*aj-T(2)*bj+T(3)*cj+AT3;
end
AT=[AT1,AT2,AT3]';
UBC=T+inv(C'*C)*AT*sgm^2;
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
Fi=Fai;
P_Bias=sqrt((UL(1)-X)^2+(UL(2)-Y)^2);
H_Bias=abs(Fi-NA)*180/pi;
OV=[X,Y,Fi];
