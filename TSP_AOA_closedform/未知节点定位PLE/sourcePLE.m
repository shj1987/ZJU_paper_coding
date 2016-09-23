clc;  %%%利用静态声源辅助实现单个未知节点定位WLS
clear all
format long
close all

%% 设定全局变量
global SIGMA
global M
%% 设定基本的参数，比如蒙特卡洛的个数，DOA误差范围，声源个数等参数
M = 5000;                         %Monte Caro times
Beacon=[10,10,pi/6;50,40,pi/5;20,60,pi/3;80,90,pi/4]; %Beacon=[10,10,pi/6;170,40,pi/5;90,110,pi/3;15,180,pi/4];
Bl=length(Beacon(:,1));           %信标个数 Beacon=[10,20,pi/3;0,160,pi/6;170,190,pi/5;130,10,pi/4;60,100,pi/5];
Unodes=[70,80,pi/3];            %未知节点的位置信息Unodes=[170,120,pi/3]; 
R=100;                            %感知半径
threshold = 1050;
err=6;
%% 在不同声源数目下进行定位
z=0;
 for SIGMA = 1: 1 : err
        num=0;
        nm=0;
        for i=1:M    
            kk=1;
            ei=0;
            pi1=0;
            SPbias=[];
 %% %%%%%%%%%%% 未知节点对声源和信标节点的DOA测量值%%%%%%%%%%%%%%%%%%%%
           S_Bcon=Beacon(:,1:2); 
            ku=0;
            for ib=1:length(S_Bcon(:,1))
                 BL=S_Bcon(ib,:);           %已定位的声源和信标
                 if sqrt((BL(1)-Unodes(1))^2+(BL(2)-Unodes(2))^2)<=R                       
                       ku=ku+1;
                       DOA(ku)= Current_Radian(BL,Unodes);   % 测得的带有噪声DOA值
                       SBStore(ku,:)=BL;    %声源估计位置存储
                 end
            end
            if ku>=3
  %% 用已知的声源对未知节点进行定位
                [OV,P_Bias, H_Bias] = PLS(SBStore,DOA,Unodes(1:2),Unodes(3));%使用Psedo-Linear Least Squares.
                [MP1, MH1]=BML1(SBStore,DOA,Unodes(1:2),Unodes(3));
                [MP2, MH2]=BML2(SBStore,DOA,Unodes(1:2),Unodes(3));
                [MP3, MH3]=BML3(SBStore,DOA,Unodes(1:2),Unodes(3),OV);
                if  P_Bias< threshold  && MP2< threshold && MP1< threshold &&  MP3< threshold 
                    num= num + 1;
                    PLS_PBias(num) = P_Bias;  %存储结果，运行次数num
                    PLS_HBias(num) = H_Bias;
                    ML1_PBias(num) = MP1;
                    ML1_HBias(num) = MH1;
                    ML2_PBias(num) = MP2;
                    ML2_HBias(num) = MH2;
                    ML3_PBias(num) = MP3;
                    ML3_HBias(num) = MH3;   
                end
            end
            DOA=[];
            SBStore=[];
        end
        display('------the programming is running now-----');
        z=z+1;
        PBIAS(z) = mean(PLS_PBias); %位置误差
        HBIAS(z)= mean(PLS_HBias);
        M1PBIAS(z) =mean(ML1_PBias); %位置误差
        M1HBIAS(z)= mean(ML1_HBias);
        M2PBIAS(z) =mean(ML2_PBias); %位置误差
        M2HBIAS(z)= mean(ML2_HBias);
        M3PBIAS(z) =mean(ML3_PBias); %位置误差
        M3HBIAS(z)= mean(ML3_HBias);
        PLS_PBias=[];
        PLS_HBias=[];
        ML1_PBias=[];
        ML1_HBias=[];
        ML2_PBias=[];
        ML2_HBias=[];
        ML3_PBias=[];
        ML3_HBias=[];
 end
MLPbias=[M1PBIAS;M2PBIAS;M3PBIAS];
MLHbias=[M1HBIAS;M2HBIAS;M3HBIAS];
Sigma=1:err;
%% 数据存储
xlswrite('PLE_PE',PBIAS);%位置误差保存
xlswrite('PLE_HE',HBIAS);%角度误差保存
xlswrite('MPBIAS',MLPbias);%位置误差保存
xlswrite('MHBIAS',MLHbias);%角度误差保存

%% 图形显示2  节点方差
figure(1)
subplot(2,1,1)
plot(Sigma,PBIAS,'b*--',Sigma, M1PBIAS,'k>--',Sigma, M2PBIAS,'rd--', Sigma, M3PBIAS,'go--','linewidth',1.5)
set(gca,'Fontsize',14);
legend('PLE','ML1','ML2','ML(ple)');
xlabel('Bearing Noise Standard Deviation (degree)');
ylabel('Location Error (m)');
xlim([1 err]);
grid on 

subplot(2,1,2)
plot(Sigma,HBIAS,'b*--', Sigma, M1HBIAS,'k>--',Sigma, M2HBIAS,'rd--', Sigma, M3HBIAS,'go--','linewidth',1.5)
set(gca,'Fontsize',14)
legend('PLE','ML1','ML2','ML(ple)');
xlabel('Bearing Noise Standard Deviation (degree)');
ylabel('Orienation Error (degree)');
xlim([1 err])
grid on 
