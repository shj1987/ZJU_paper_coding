function [EL,EH]=BMLThy_Node(TSB,AOA,Unode)
global SIGMA
Sgm=SIGMA*pi/180;
W=length(AOA);
F11=0;
F12=0;
F13=0;
F21=0;
F22=0;
F23=0;
F31=0;
F32=0;
F33=0;
for ia=1:W
    HP=TSB(ia,:);
    Dr=(Unode(1)-HP(1))^2+(Unode(2)-HP(2))^2;
    Sgd=1/(Dr^2*Sgm^2);
    F11=F11+Sgd*(Unode(2)-HP(2))^2;
    F12=F12-Sgd*(Unode(1)-HP(1))*(Unode(2)-HP(2));
    F13=F13+Sgd*Dr*(Unode(2)-HP(2));
    F21=F21-Sgd*(Unode(1)-HP(1))*(Unode(2)-HP(2));
    F22=F22+(Unode(1)-HP(1))^2;
    F23=F23-Sgd*Dr*(Unode(1)-HP(1));
    F31=F31+Sgd*Dr*(Unode(2)-HP(2));
    F32=F32-Sgd*Dr*(Unode(1)-HP(1));
    F33=F33+Sgm^2;
end
F=[F11 F12 F13;F21 F22 F23;F31 F32 F33];
FI=inv(F);
EL=sqrt(FI(1,1)+FI(2,2));
EH=sqrt(FI(3,3));



