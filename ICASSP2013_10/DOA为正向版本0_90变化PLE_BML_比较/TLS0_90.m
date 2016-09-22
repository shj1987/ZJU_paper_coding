clc;  %%%利用静态声源辅助实现单个未知节点定位WLS
clear all
format long
close all

%% 设定全局变量
global SIGMA
global M
%% 设定基本的参数，比如蒙特卡洛的个数，DOA误差范围，声源个数等参数
M = 2000;                         %Monte Caro times
R=200;                            %感知半径
threshold = 200;
err=180;
SIGMA=4;
%% 在不同声源数目下进行定位
z=0;
 for Hd = 0: 2: err
        Hede=Hd*pi/180;
        %Unodes=[Unodep,Hede];
        num=0;
        nm=0;
        for i=1:M    
            kk=1;
            ei=0;
            pi1=0;
            SPbias=[];
            Beacon=round(rand(8,2)*100);
            node=round(rand(1,2)*100);
%             Hdr=round(rand(1,1)*100*pi/2)/100;
            Unodes=[node,Hede];
 %% %%%%%%%%%%% 未知节点对声源和信标节点的DOA测量值%%%%%%%%%%%%%%%%%%%%
           S_Bcon=Beacon; 
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
                [OV,PLE_PBias, PLE_HBias] = PLS(SBStore,DOA,Unodes(1:2),Unodes(3)); %PLE
                [OV2,P_Bias, H_Bias] = AVTLS(SBStore,DOA,Unodes(1:2),Unodes(3));%使用BCPLE
%                 [WP_Bias, WH_Bias] = WIV(SBStore,DOA,Unodes(1:2),Unodes(3),U);%%%BCAVPLEWIV
%                 [PLE_WP, PLE_WH] = WIV2(SBStore,DOA,Unodes(1:2),Unodes(3),UT);%AVPLE_WIV
                [TP_Bias, TH_Bias] = TRI_TWO(SBStore,DOA,Unodes(1:2),Unodes(3));
                [MP, MH]=BMLF(SBStore,DOA,Unodes(1:2),Unodes(3),OV2);
                if  PLE_PBias< threshold && P_Bias< threshold && MP<threshold && TP_Bias<threshold
                    num= num + 1;
                    TP_B(num)=TP_Bias;
                    TH_B(num)=TH_Bias;
                    PE_Pbias(num)=PLE_PBias;
                    PE_Hbias(num)=PLE_HBias;
                    BC_PBias(num) = P_Bias;  %存储结果，运行次数num
                    BC_HBias(num) = H_Bias;
                    ML_PBias(num) = MP;
                    ML_HBias(num) = MH; 
                end
            end
            DOA=[];
            SBStore=[];
        end
        display('------the programming is running now-----');
        z=z+1;
        TP(z)=mean(TP_B);
        TH(z)=mean(TH_B);
        PLE_P(z)=mean(PE_Pbias);
        PLE_H(z)=mean(PE_Hbias);
        PBIAS(z) = mean(BC_PBias); %位置误差
        HBIAS(z)= mean(BC_HBias);
        MPBIAS(z) =mean(ML_PBias); %位置误差
        MHBIAS(z)= mean(ML_HBias);
        BC_PBias=[];
        BC_HBias=[];
        ML_PBias=[];
        ML_HBias=[];
        WIV_Pbias=[];
        WIV_Hbias=[];
        PE_Pbias=[];
        PE_Hbias=[];
        PLE_WPbias=[];
        PLE_WHbias=[];
 end
Sigma = 0: 2 : err;
PB=[PLE_P;PBIAS;MPBIAS];%ple, bcple,ML. avplewiv, bcplvwiv
HB=[PLE_H;HBIAS;MHBIAS];
xlswrite('PB_90dg', PB);        %位置误差保存
xlswrite('PH_90dg', HB);        %角度误差保存

%% 图形显示2  节点方差
figure(1)
subplot(2,1,1)
plot(Sigma,TP,'kx--',Sigma,PLE_P,'m+-',Sigma,PBIAS,'b.--', Sigma, MPBIAS,'g.--','linewidth',1.5)
%plot(Sigma,PLE_P,'m--',Sigma,PBIAS,'b-.', Sigma, MPBIAS,'r*--', Sigma, WIV_PB,'g+--','linewidth',1.5)
set(gca,'Fontsize',13);
h=legend('Triangulation','AVPLE','AVTLS','ML');
set(h,'Fontsize',12)
xlabel('Orientation of Unknown Node (degrees)');
ylabel('Location Error (m)');
axis([0 180 0 40]);
set(gca,'XTick',0:20:180);
grid on 

subplot(2,1,2)
plot(Sigma,TH,'kx--',Sigma,PLE_H,'m+-',Sigma,HBIAS,'b.--',Sigma, MHBIAS,'g.--','linewidth',1.5)
%plot(Sigma,PLE_H,'m>--',Sigma,HBIAS,'b*--', Sigma, WIV_HB,'g+--',Sigma, PLE_WIVH,'k^--','linewidth',1.5)
set(gca,'Fontsize',13)
h=legend('Triangulation','AVPLE','AVTLS','ML');
set(h,'Fontsize',12)
xlabel('Orientation of Unknown Node (degrees)');
ylabel('Orienation Error (degrees)');
axis([0 180 0 50]);
set(gca,'XTick',0:20:180);
grid on 
