clc;  %%%利用静态声源辅助实现单个未知节点定位
clear all
format long
close all

%% 设定全局变量
display('----start now-----');
global SIGMA
global M
%% 设定基本的参数，比如蒙特卡洛的个数，DOA误差范围，声源个数等参数
M = 10000;                %Monte Caro times
NumSource =50;            %unknown node的个数
NumB=30;                  %信标节点
R=100;                    %感知半径
threshold = 150;
err=4;
%% 在不同声源数目下进行定位
 for ii=1:6
    ASource=NumSource;
    z=0;
    for SIGMA = 3: 1 : 3
        nu=0;
        nm=0;
        for i=1:M           
            SL1=round(rand(ASource/2,2)*500); %多个声源位置分布500*500
            SL2=round(rand(ASource/2,2)*500);
            SL=[SL1;SL2];
            UA=round(rand(ASource,1)*pi*10)/10;  
            SLocation=[SL,UA];
            BL1=round(rand(NumB/2,2)*500);
            BL2=round(rand(NumB/2,2)*500);
            BL=[BL1;BL2];
            BA=round(rand(NumB,1)*pi*10)/10;
            Beacon=[BL,BA];
            kk=1;
            ei=0;
            pi1=0;
            SPbias=[];
            S_Source=[];
            A_Source=[];
            jj=0;
            Unode=[];
            cnt=0;
     %% 信标定位静态声源与未知节点%%%%%%%%
            for j=1: ASource
                for k=1:NumB
                    BBeacon=Beacon(k,:);           %信标信息
                    SSource=SLocation(j,:);  %声源信息
                    if sqrt((SSource(1)-BBeacon(1))^2+(SSource(2)-BBeacon(2))^2)<=R   %感知半径范围之内
                       XITA= DOAMeasure(BBeacon,SSource); % 测声源的DOA值子函数
                       XXITA(kk)=XITA;         %绝对坐标系下的角度
                       Beacon_use(kk,:)= BBeacon;          %能够测到的信标节点
                       kk=kk+1;
                    end
                end
     %% 信标节点对静态声源进行定位
                 if (kk-1)>=3
                     nu=nu+1;
                     Deg(nu)=kk-1;
                    [SPT,pc]=PLS1(Beacon_use(:,1:2),XXITA,SSource(1:2),SSource(3));
                    if pc(1)<threshold  %设定阈值
                       ei=ei+1;
                       pi1=pi1+1;
                       SPbias(pi1)=pc(1);
                       S_Source(ei,:)=SPT;  %已经定位的声源存储
                       A_Source(ei,:)=SSource(1:2);%%准确位置
                       cnt=cnt+1;
                    else jj=jj+1; Unode(jj,:)=SSource; 
                    end
                 else
                     jj=jj+1;
                     Unode(jj,:)=SSource;     
                 end
                 XXITA=[];
                 Beacon_use=[];   
                 kk=1;             
            end
 %% %%%%%%%%%%% 未知节点对声源和信标节点的DOA测量值%%%%%%%%%%%%%%%%%%%%
            S_Bcon=Beacon(:,1:2);
            B_Source=[];
            AB_Source=[];
            B_Source=[S_Source;S_Bcon];         %信标和声源位置信息
            AB_Source=[A_Source;S_Bcon];      %准确位置  
            if isempty(Unode)==1
                NumU=0;
            else
               NumU=length(Unode(:,1));
            end
            for ij=1:NumU
                ku=0;
                for iu=1:length(B_Source(:,1))
                     Unodes=[];
                     BSL=AB_Source(iu,:);           %已定位的声源和信标
                     E_BSL=B_Source(iu,:);          %估计位置
                     Unodes=Unode(ij,:);
                     if sqrt((BSL(1)-Unodes(1))^2+(BSL(2)-Unodes(2))^2)<=R                       
                           ku=ku+1;
                           DOA(ku)= DOAMeasure(BSL,Unodes);   % 测得的带有噪声DOA值
                           SBStore(ku,:)=E_BSL;    %声源估计位置存储
                     end
                end
                nu=nu+1;
                Deg(nu)=ku;
                DOA=[];
                SBStore=[];
                S_Source=[];  %一次结束
            end
        end
        display('---------the programing is running now-----------------');
        z=z+1;
        Num(ii,z)=mean(Deg); %参与定位个数
        SPavg=[];
        Deg=[];
   end
    
 end
%% 数据存储
xlswrite('NodeD30',Num);%位置误差保存
display('-----OK***********NOW----')


