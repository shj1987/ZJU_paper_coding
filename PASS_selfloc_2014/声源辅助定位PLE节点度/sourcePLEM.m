clc;  %%%利用静态声源辅助实现单个未知节点定位
clear all
format long
close all

%% 设定全局变量
global SIGMA
global M
%% 设定基本的参数，比如蒙特卡洛的个数，DOA误差范围，声源个数等参数
M = 1000;                         %Monte Caro times
NumSource =[40 60 80 100 120];            %声源的个数
NumB=30;                           %信标节点
NumU=60;                           %未知节点个数 
R=100;                             %感知半径
threshold = 100;
err=3;
z=0;
%% 在不同声源数目下进行定位
 for ii=1:length(NumSource)
    ASource=NumSource(ii);
    for SIGMA = 2: 1 : 2
        nm=0;
        for i=1:M           
            SL1=round(rand(ASource,2)*500); %多个声源位置分布500*500
            %SL2=round(rand(ASource/2,2)*500);
            SLocation=SL1;
            BL=round(rand(NumB,2)*500);
%             BL2=round(rand(NumB/2,2)*500);
%             BL=[BL1;BL2];
            BA=round(rand(NumB,1)*pi*10)/10;
            Beacon=[BL,BA];
            UL=round(rand(NumU,2)*500);        %未知节点
%             UL2=round(rand(NumU/2,2)*500);        %未知节点
%             UL=[UL1;UL2];
            UA=round(rand(NumU,1)*pi*10)/10;  
            Unode=[UL,UA];
            kk=1;
            ei=0;
            pi1=0;
            SPbias=[];
            S_Source=[];
            A_Source=[];
     %% 信标定位静态声源与未知节点%%%%%%%%
            for j=1: ASource
                for k=1:NumB
                    BBeacon=Beacon(k,:);           %信标信息
                    SSource=SLocation(j,:);        %声源信息
                    if sqrt((SSource(1)-BBeacon(1))^2+(SSource(2)-BBeacon(2))^2)<=R   %感知半径范围之内
                       XITA= DOAMeasure(SSource,BBeacon); % 测声源的DOA值子函数
                       XXITA(kk)=XITA+BBeacon(3);             %绝对坐标系下的角度
                       Beacon_use(kk,:)= BBeacon;             %能够测到的信标节点
                       kk=kk+1;
                    end
                end
     %% 信标节点对静态声源进行定位
                 if (kk-1)>=2
                    OV1 = TLS_Single(XXITA,Beacon_use(:,1:2));    %声源估计位置存储PLE
                    Spt = BML1(XXITA,Beacon_use(:,1:2),OV1);
                    pc=sqrt((SSource(1)-Spt(1))^2+(SSource(2)-Spt(2))^2); %位置偏差
                    if all(OV1==0)~=1 && pc<threshold  %设定阈值
                       ei=ei+1;
                       pi1=pi1+1;
                       SPbias(pi1)=pc;
                       S_Source(ei,:)=Spt;  %已经定位的声源存储
                       A_Source(ei,:)=SSource;%%准确位置
                    end
                 end
                 XXITA=[];
                 Beacon_use=[];   
                 kk=1;             
            end
            
 %% %%%%%%%%%%% 未知节点对声源和信标节点的DOA测量值%%%%%%%%%%%%%%%%%%%%
            S_Bcon=Beacon(:,1:2);
            B_Source=[];
            AB_Source=[];
            B_Source=[S_Source;S_Bcon];       %信标和声源位置信息
            AB_Source=[A_Source;S_Bcon];      %准确位置  
            cnt=0;
            nu1=0;
            nu2=0;
            nu3=0;
            nu4=0;
            nu5=0;
            for ij=1:NumU
                ku=0;
                for iu=1:length(B_Source(:,1))
                     %Unodes=[];
                     BSL=AB_Source(iu,:);           %已定位的声源和信标
                     E_BSL=B_Source(iu,:);          %估计位置
                     Unodes=Unode(ij,:);
                     if sqrt((BSL(1)-Unodes(1))^2+(BSL(2)-Unodes(2))^2)<=R                       
                           ku=ku+1;
                           DOA(ku)= DOAMeasure(BSL,Unodes);   % 测得的带有噪声DOA值
                           SBStore(ku,:)=E_BSL;    %声源估计位置存储
                     end
                end
%                 if ku>=3
%                     nu1=nu1+1;
%                 end
%                 if ku>=4
%                     nu2=nu2+1;
%                 end
%                 if ku>=5
%                     nu3=nu3+1;
%                 end
%                 if ku>=6
%                     nu4=nu4+1;
%                 end
%                 if ku>=7
%                     nu5=nu5+1;
%                 end
                DOA=[];
                SBStore=[];
                S_Source=[];  %一次结束
            end
            SKU(i)=ku;
        end       
        display('---------the programing is running now-----------------');
        z=z+1;
        NKU(z)=mean(SKU);
        SPavg=[];
        SKU=[];
        PI=[];
   end
    
 end
%% 数据存储
xlswrite('NodeDegree',NKU);%位置误差保存
display('----******-OK--*******--');
NKU


