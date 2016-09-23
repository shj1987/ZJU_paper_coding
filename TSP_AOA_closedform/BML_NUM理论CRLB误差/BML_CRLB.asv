clc;  %%%利用静态声源辅助实现单个未知节点定位WLS
clear all
format long
close all
%%%%% --------CRLB-------- 误差%%%%
%% 设定全局变量
global SIGMA
global M
%% 设定基本的参数，比如蒙特卡洛的个数，DOA误差范围，声源个数等参数
M = 2000;                         %Monte Caro times
%Beacon=[45,55;30,70;95,50;100,50;85,35;40,40;50,15;60,40;20,30;75,60];
% Beacon=round(rand(100,2)*100);
% Unodes=[100,90,pi/3];            %未知节点的位置信息Unodes=[170,120,pi/3]; 
R=200;                            %感知半径
threshold = 150;
err=3.5;
%% 在不同声源数目下进行定位
z=0;
nm=0;
SIGMA=2;
 for i=1:M 
        num=0;
        nm=nm+1;
        node=round(rand(1,2)*100);
        Hdr=round(rand(1,1)*10*pi/2)/10;
        Unodes=[node,Hdr];
        for NUM = 10: 50: 410  
            kk=1;
            ei=0;
            pi1=0;
            SPbias=[];
            Beacon=round(rand(NUM,2)*100);
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
                [U,OV,BP_Bias, BH_Bias] = BCPLS(SBStore,DOA,Unodes(1:2),Unodes(3));%使用Psedo-Linear Least Squares.
                [MP, MH]=BMLF(SBStore,DOA,Unodes(1:2),Unodes(3),OV);
                if  MP< threshold 
                    [TP_Bias,TH_Bias]=BMLThy_Node(SBStore,DOA,Unodes);
                    num= num + 1;
                    BCPLS_PBias(nm,num) = BP_Bias;  %存储结果，运行次数num
                    BCPLS_HBias(nm,num) = BH_Bias;
                    ML_PBias(nm,num) = MP;
                    ML_HBias(nm,num) = MH; 
                    CRLB_Pbias(nm,num)= TP_Bias;
                    CRLB_Hbias(nm,num)= TH_Bias;
                else nm=nm-1;
                    break;
                end
            end
            DOA=[];
            SBStore=[];
        end
        display('------the programming is running now-----');
%         z=z+1;
%         PBIAS(z) = mean(BCPLS_PBias); %位置误差
%         HBIAS(z)= mean(BCPLS_HBias);
%         MPBIAS(z) =mean(ML_PBias); %位置误差
%         MHBIAS(z)= mean(ML_HBias);
%         CRLB_PB(z)=mean(CRLB_Pbias);
%         CRLB_HB(z)=mean(CRLB_Hbias);
%         BCPLS_PBias=[];
%         BCPLS_HBias=[];
%         ML_PBias=[];
%         ML_HBias=[];
%         WIV_Pbias=[];
%         WIV_Hbias=[];
 end
PBIAS = mean(BCPLS_PBias); %位置误差
HBIAS= mean(BCPLS_HBias);
MPBIAS=mean(ML_PBias); %位置误差
MHBIAS= mean(ML_HBias);
CRLB_PB=mean(CRLB_Pbias);
CRLB_HB=mean(CRLB_Hbias);
 NUM = 10: 50: 410 ;
%% 数据存储
xlswrite('CRLB_PB410',CRLB_PB);%ML位置误差保存
xlswrite('CRLB_HB410',CRLB_HB);%ML角度误差保存


%% 图形显示2  节点方差
figure(1)
subplot(2,1,1)
plot(NUM,CRLB_PB,'b*--','linewidth',1.5)
set(gca,'Fontsize',14);
% legend('BCPLE','ML','CRLB');
xlabel('Bearing Noise Standard Deviation (degree)');
ylabel('Location Error (m)');
grid on 

subplot(2,1,2)
plot(NUM,CRLB_PB,'b*--','linewidth',1.5)
set(gca,'Fontsize',14)
% legend('BCPLE','ML','CRLB');
xlabel('Bearing Noise Standard Deviation (degree)');
ylabel('Orienation Error (degree)');
grid on 
