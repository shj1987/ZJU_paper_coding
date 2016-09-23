clc;  %%%利用静态声源辅助实现单个未知节点定位WLS
clear all
format long
close all

%% 设定全局变量
global SIGMA
global M
%% 设定基本的参数，比如蒙特卡洛的个数，DOA误差范围，声源个数等参数
M = 2000;                         %Monte Caro times
%Beacon=[80,55;30,70;95,50;100,50;85,35;0,40;50,15;60,40;20,30;75,60];
%Beacon=[45,55;30,70;95,50;100,50;85,35;40,40;50,15;60,40;20,30;75,60];
%Bl=length(Beacon(:,1));           %信标个数
%Unodes=[100,90,pi/3];
% plot(Beacon(:,1),Beacon(:,2),'r>');
% plot(Unodes(1),Unodes(2));
R=200;                            %感知半径
threshold = 150;
err=10;
%% 在不同声源数目下进行定位
z=0;
SIGMA=2;
 for BN=5: 5: 40
        num=0;
        nm=0;
        for i=1:M    
            kk=1;
            ei=0;
            pi1=0;
            SPbias=[];
            Beacon=round(rand(BN,2)*100);
            node=round(rand(1,2)*100);
            Hdr=round(rand(1,1)*100*pi/2)/100;
            Unodes=[node,Hdr];
 %% %%%%%%%%%%% 未知节点对声源和信标节点的DOA测量值%%%%%%%%%%%%%%%%%%%%
           S_Bcon=Beacon(:,1:2); 
            ku=0;
            for ib=1:length(S_Bcon(:,1))
                 BL=S_Bcon(ib,:);           %已定位的声源和信标
                 if sqrt((BL(1)-Unodes(1))^2+(BL(2)-Unodes(2))^2)<=R                       
                       ku=ku+1;
                       DOA1(ku)= Current_Radian(BL,Unodes);   % 测得的带有噪声DOA值
                       DOA2(ku)= Current_Radian2(BL,Unodes);
                       SBStore(ku,:)=BL;    %声源估计位置存储
                 end
            end
            if ku>=3
  %% 用已知的声源对未知节点进行定位
                [UT1,PLE_PBias, PLE_HBias] = PLS(SBStore,DOA1,Unodes(1:2),Unodes(3)); %PLE
                [UT2,PLE_PBias, PLE_HBias] = PLS(SBStore,DOA2,Unodes(1:2),Unodes(3)); %PLE
%                 [U,OV,P_Bias, H_Bias] = BCPLS(SBStore,DOA,Unodes(1:2),Unodes(3));%使用BCPLE
%                 [WP_Bias, WH_Bias] = WIV(SBStore,DOA,Unodes(1:2),Unodes(3),U);%%%BCAVPLEWIV
                [PLE_WP, PLE_WH] = WIV2(SBStore,DOA1,Unodes(1:2),Unodes(3),UT1);%AVPLE_WIV
                [PLE_WP2, PLE_WH2] = PLEWIV(SBStore,DOA2,Unodes(1:2),Unodes(3),UT2);
                %[MP, MH]=BMLF(SBStore,DOA,Unodes(1:2),Unodes(3),OV);
                if  PLE_WP2< threshold && PLE_WP< threshold
                    num= num + 1;
                    PE_Pbias(num)=PLE_PBias;
                    PE_Hbias(num)=PLE_HBias;
%                     BC_PBias(num) = P_Bias;  %存储结果，运行次数num
%                     BC_HBias(num) = H_Bias;
                    %ML_PBias(num) = MP;
                    %ML_HBias(num) = MH; 
                    PLE_WPbias(num)=PLE_WP;
                    PLE_WHbias(num)=PLE_WH;
                    WIV_Pbias(num)= PLE_WP2;
                    WIV_Hbias(num)= PLE_WH2;
                end
            end
            DOA1=[];
            DOA2=[];
            SBStore=[];
        end
        display('------the programming is running now-----');
        z=z+1;
        PLE_P(z)=mean(PE_Pbias);
        PLE_H(z)=mean(PE_Hbias);
        WIV_PB(z)=mean(WIV_Pbias);
        WIV_HB(z)=mean(WIV_Hbias);
        PLE_WIVP(z)=mean(PLE_WPbias);
        PLE_WIVH(z)=mean(PLE_WHbias);
        WIV_Pbias=[];
        WIV_Hbias=[];
        PE_Pbias=[];
        PE_Hbias=[];
        PLE_WPbias=[];
        PLE_WHbias=[];
 end
Sigma = 5: 5: 40;
%% 数据存储
PB=[PLE_WIVP;WIV_PB];%ple, bcple, ML, avplewiv, bcplvwiv
HB=[PLE_WIVH;WIV_HB];
xlswrite('PBwiv', PB);        %位置误差保存
xlswrite('PHwiv', HB);        %角度误差保存

% xlswrite('PLE_PE10', PLE_P);        %位置误差保存
% xlswrite('PLE_HE10', PLE_H);        %角度误差保存
% xlswrite('BCPLE_PE10',PBIAS);       %位置误差保存
% xlswrite('BCPLE_HE10',HBIAS);       %角度误差保存
% xlswrite('MLPB10',MPBIAS);          %ML位置误差保存
% xlswrite('MLHB10',MHBIAS);          %ML角度误差保存
% xlswrite('BCPLEWIV_PB10',WIV_PB);   %ML位置误差保存
% xlswrite('BCPLEWIV_HB10',WIV_HB);   %ML角度误差保存
% xlswrite('PLEWIV_PB10', PLE_WIVP);   %ML位置误差保存
% xlswrite('PLEWIV_HB10', PLE_WIVH);   %ML角度误差保存


%% 图形显示2  节点方差
figure(1)
subplot(2,1,1)
%plot(Sigma,PLE_P,'m>--',Sigma,PBIAS,'b*--', Sigma, MPBIAS,'rs--', Sigma, WIV_PB,'g+--',Sigma, PLE_WIVP,'ko--','linewidth',1.5)
%plot(Sigma,PLE_P,'m--',Sigma,PBIAS,'b-.', Sigma, MPBIAS,'r*--', Sigma, WIV_PB,'g+--','linewidth',1.5)
plot(Sigma,WIV_PB,'k<--',Sigma, PLE_WIVP,'ro--','linewidth',1.5)
set(gca,'Fontsize',13);
h=legend('AVPLE-WIV(not i.i.d)','AVPLE-WIV (i.i.d)');
set(h,'Fontsize',12)
xlabel('Number of Beacons (N)');
ylabel('Location Error (m)');
% set(gca,'XTick',1:10);
% axis([1 10 0 20])
grid on 

subplot(2,1,2)
%plot(Sigma,PLE_H,'m>--',Sigma,HBIAS,'b*--', Sigma, MHBIAS,'rs--', Sigma, WIV_HB,'g+--',Sigma, PLE_WIVH,'ko--','linewidth',1.5)
%plot(Sigma,PLE_H,'m>--',Sigma,HBIAS,'b*--', Sigma, WIV_HB,'g+--',Sigma, PLE_WIVH,'k^--','linewidth',1.5)
plot(Sigma, WIV_HB,'k<--',Sigma, PLE_WIVH,'ro--','linewidth',1.5)
set(gca,'Fontsize',13)
%legend('PLE','BCPLE','ML','BCAVPLE-WIV','AVPLE-WIV');
h=legend('AVPLE-WIV(not i.i.d)','AVPLE-WIV (i.i.d)');
set(h,'Fontsize',12)
xlabel('Number of Beacons (N)');
ylabel('Orienation Error (degrees)');
% set(gca,'XTick',1:10);
% axis([1 10 0 30])
grid on 
