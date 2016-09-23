clc;  %%%利用静态声源辅助实现单个未知节点定位BCPLE
clear all
format long
close all
%% 设定全局变量
global SIGMA
%% 设定基本的参数，比如蒙特卡洛的个数，DOA误差范围，声源个数等参数
M = 2000;                         %Monte Caro times
%Beacon=[50 60,pi/10;60,30,pi/6;10,10,pi/6;50,40,pi/5;20,60,pi/3;80,90,pi/4]; %Beacon=[10,10,pi/6;170,40,pi/5;90,110,pi/3;15,180,pi/4];
Beacon=rand(10,2)*100;
Bl=length(Beacon(:,1));           %信标个数 Beacon=[10,20,pi/3;0,160,pi/6;170,190,pi/5;130,10,pi/4;60,100,pi/5];
Unodes=[50,80,pi/3];            %未知节点的位置信息Unodes=[170,120,pi/3]; 
R=300;                            %感知半径
threshold = 200;
err=6;
%% 在不同声源数目下进行定位
z=0;
 for SIGMA = 0.5: 0.5 : 4
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
                [P_Bias, H_Bias] = PLS(SBStore,DOA,Unodes(1:2),Unodes(3));%使用Psedo-Linear Least Squares.
                [B_xyz,BP_Bias, BH_Bias] = BCPLS(SBStore,DOA,Unodes(1:2),Unodes(3));%使用Bias Compensated PseudoLinear Least Squares.
                if  P_Bias< threshold
                    num= num + 1;
                    PLS_PBias(num) = P_Bias;  %存储结果，运行次数num
                    PLS_HBias(num) = H_Bias;
                    BCPL_PBias(num)=BP_Bias;
                    BCPL_HBias(num)=BH_Bias;
                end
            end
            DOA=[];
            SBStore=[];
        end
        display('------the programming is running now-----');
        z=z+1;
        PBIAS(z)= mean(PLS_PBias); %位置误差
        HBIAS(z)= mean(PLS_HBias);
        BC_PBIAS(z)= mean(BCPL_PBias);
        BC_HBIAS(z)= mean(BCPL_HBias);
        PLS_PBias=[];
        PLS_HBias=[];
        BCPL_PBias=[];
        BCPL_HBias=[];
 end

Sigma=0.5: 0.5 : 4;
%% 数据存储
xlswrite('BCPLE_PE',PBIAS);  %位置误差保存
xlswrite('BCPLE_HE',HBIAS);  %角度误差保存
xlswrite('BCPLE_PE',BC_PBIAS); %位置误差保存
xlswrite('BCPLE_HE',BC_HBIAS); %角度误差保存

%% 图形显示2  节点方差
figure(1)
subplot(2,1,1)
plot(Sigma, PBIAS,'b*--',Sigma, BC_PBIAS,'go--','linewidth',1.5)
set(gca,'Fontsize',14);
legend('PLE','BCPLE');
xlabel('Bearing Noise Standard Deviation (degree)');
ylabel('Location Error (m)');
xlim([0.5 4]);
grid on 

subplot(2,1,2)
plot(Sigma, HBIAS,'b*--',Sigma, BC_HBIAS,'go--','linewidth',1.5)
set(gca,'Fontsize',14)
legend('PLE','BCPLE');
xlabel('Bearing Noise Standard Deviation (degree)');
ylabel('Orienation Error (degree)');
xlim([0.5 4]);
grid on 
