function DOAS=DOA2_Measure(L,SNode,Node) %beaacon, neighboring node and unknown node
global SIGMA
A=Node;
B=SNode(1,1:2);
CC=SNode(2,1:2);
XT1=atan2(B(2)-Node(2), B(1)-Node(1))+normrnd(0, SIGMA)*pi/180-Node(3);
XT2=atan2(CC(2)-Node(2), CC(1)-Node(1))+normrnd(0, SIGMA)*pi/180-Node(3);
if XT1<XT2
    B=SNode(2,1:2);
    CC=SNode(1,1:2);
end
BL = atan2(L(2)-B(2), L(1)-B(1))+normrnd(0, SIGMA)*pi/180;    %绝对坐标系的值在[-180，180]之间根据象限确定
if BL<0
   BL=BL+2*pi;
end
BC = atan2(CC(2)-B(2), CC(1)-B(1))+normrnd(0, SIGMA)*pi/180;  
if BC<0
   BC=BC+2*pi;
end
CBL =BL-BC;  %angle CBL
if CBL<0
   CBL=CBL+2*pi;
end
%%%%%%%%%%%%%%%
CB = atan2(B(2)-CC(2), B(1)-CC(1))+normrnd(0, SIGMA)*pi/180;
if CB<0
   CB=CB+2*pi;
end
CL= atan2(L(2)-CC(2), L(1)-CC(1))+normrnd(0, SIGMA)*pi/180;
if CL<0
    CL=CL+2*pi;
end
BCL=CB-CL;
if BCL<0
    BCL=BCL+2*pi;
end
%%%%%%%%% ABC %%%%%%%
BA = atan2(A(2)-B(2), A(1)-B(1))+normrnd(0, SIGMA)*pi/180;
if BA<0
    BA=BA+2*pi;
end
ABC =BC-BA;
if ABC<0
    ABC=ABC+2*pi;
end
%%%%% BAC %%%%%%
AB=  atan2(B(2)-A(2), B(1)-A(1))+normrnd(0, SIGMA)*pi/180;
if AB<0
    AB=AB+2*pi;
end
AC = atan2(CC(2)-A(2), CC(1)-A(1))+normrnd(0, SIGMA)*pi/180;
if AC<0
    AC=AC+2*pi;
end
BAC=AB-AC;
if BAC<0
    BAC=BAC+2*pi;
end
%%% ACB %%%%
CA = atan2(A(2)-CC(2), A(1)-CC(1))+normrnd(0, SIGMA)*pi/180;
if CA<0
    CA=CA+2*pi;
end
ACB = CA-CB;
if ACB<0
    ACB=ACB+2*pi;
end
%%% BLC%%%%
LB= atan2(B(2)-L(2), B(1)-L(1))+normrnd(0, SIGMA)*pi/180;
if LB<0
    LB=LB+2*pi;
end
LC= atan2(CC(2)-L(2), CC(1)-L(1))+normrnd(0, SIGMA)*pi/180;
if LC<0
    LC=LC+2*pi;
end
BLC=LC-LB;
if BLC<0
    BLC=BLC+2*pi;
end
%%% ACL%%%%%
ACL=CA-CL;
if ACL<0
    ACL=ACL+2*pi;
end
%%% length %%%%%%
DCL=sin(CBL)/sin(BCL);
DBC=sin(BLC)/sin(BCL);
DAC=DBC*sin(ABC)/sin(BAC);
DAL=sqrt(DAC^2+DCL^2-2*DAC*DCL*cos(ACL));
SLAC=DCL*sin(ACL)/DAL;
LAC=asin(SLAC);
OAC=AC-Node(3);  %%DOA
if OAC<0
    OAC=OAC+2*pi;
end
DOATR= atan2(L(2)-Node(2), L(1)-Node(1))-Node(3);  %true angle
if DOATR<0
   DOATR=DOATR+2*pi;
end
DOAS=OAC+LAC;
if abs(DOATR-DOAS)>2
    DOAS=[];
end

