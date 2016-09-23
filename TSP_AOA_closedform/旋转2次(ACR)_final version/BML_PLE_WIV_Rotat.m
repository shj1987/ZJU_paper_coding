clc;  %%%利用静态声源辅助实现单个未知节点定位WLS
clear all
format long
close all

%% 设定全局变量
global SIGMA
global M
%% 设定基本的参数，比如蒙特卡洛的个数，DOA误差范围，声源个数等参数
M = 2000;                         %Monte Caro times
%Beacon=[45,55;30,70;95,50;100,50;85,35;40,40;50,15;60,40;20,30;75,60];
%Bl=length(Beacon(:,1));           %信标个数
%node=[80,25];
% plot(Beacon(:,1),Beacon(:,2),'r>');
% plot(Unodes(1),Unodes(2));
R=200;                            %感知半径
threshold = 150;
err=3.5;
SIGMA=2;
%% 在不同声源数目下进行定位
z=0;
 for Hd =0: 2: 180
        Hede=Hd*pi/180;
        %Unodes=[Unodep,Hede];
        num=0;
        nm=0;
        for i=1:M    
            kk=1;
            ei=0;
            pi1=0;
            SPbias=[];
            Beacon=round(rand(10,2)*100);
            node=round(rand(1,2)*100);
%           Hdr=round(rand(1,1)*100*pi/2)/100;
            Unodes=[node,Hede];
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
                [UT,RD1,PLE_PBias1, PLE_HBias1] = PLS(SBStore,DOA,Unodes(1:2),Unodes(3)); %PLE
                [RT,RD2,PLE_PBias2,PLE_HBias2,zg] = RPLS(SBStore,DOA,Unodes(1:2),Unodes(3));
                k1=0; 
                k2=0;
%                 for kj=1:length(RD1)
%                     if abs(RD1(kj))>abs(RD2(kj)) %error
%                         k1=k1+1;
%                     else
%                         k2=k2+1;
%                     end
%                 end
                for kj=1:length(RD1)
                   k1=(RD1(kj))^2+k1;
                   k2=(RD2(kj))^2+k2;
                end
                flg=0;
                if k1>k2 %&& zg==0;% not 180 degree
                    fg=1;
                else
                    fg=0;
                end
                if fg==1
                    PLE_PBias=PLE_PBias2;
                    PLE_HBias=PLE_HBias2;
                    [RU,OV2,P_Bias, H_Bias] = RBCPLS(SBStore,DOA,Unodes(1:2),Unodes(3));%使用BCPLE
                    [WP_Bias, WH_Bias] = RWIV(SBStore,DOA,Unodes(1:2),Unodes(3),RU);%%%BCAVPLEWIV
                    [PLE_WP, PLE_WH] = RWIV2(SBStore,DOA,Unodes(1:2),Unodes(3),RT);%AVPLE_WIV
                    [MP, MH]=BMLF(SBStore,DOA,Unodes(1:2),Unodes(3),OV2);
                else
                    PLE_PBias=PLE_PBias1;
                    PLE_HBias=PLE_HBias1;
                    [U,OV,P_Bias, H_Bias] = BCPLS(SBStore,DOA,Unodes(1:2),Unodes(3));%使用BCPLE
                    [WP_Bias, WH_Bias] = WIV(SBStore,DOA,Unodes(1:2),Unodes(3),U);%%%BCAVPLEWIV
                    [PLE_WP, PLE_WH] = WIV2(SBStore,DOA,Unodes(1:2),Unodes(3),UT);%AVPLE_WIV
                    [MP, MH]=BMLF(SBStore,DOA,Unodes(1:2),Unodes(3),OV);
                end
                if  PLE_PBias< threshold && MP<100 && MH<100
                    num= num + 1;
                    PE_Pbias(num)=PLE_PBias;
                    PE_Hbias(num)=PLE_HBias;
                    BC_PBias(num) = P_Bias;  %存储结果，运行次数num
                    BC_HBias(num) = H_Bias;
                    ML_PBias(num) = MP;
                    ML_HBias(num) = MH; 
                    PLE_WPbias(num)=PLE_WP;
                    PLE_WHbias(num)=PLE_WH;
                    WIV_Pbias(num)= WP_Bias;
                    WIV_Hbias(num)= WH_Bias;
                end
            end
            DOA=[];
            SBStore=[];
        end
        display('------the programming is running now-----');
        z=z+1;
        PLE_P(z)=mean(PE_Pbias);
        PLE_H(z)=mean(PE_Hbias);
        PBIAS(z) = mean(BC_PBias); %位置误差
        HBIAS(z)= mean(BC_HBias);
        MPBIAS(z) =mean(ML_PBias); %位置误差
        MHBIAS(z)= mean(ML_HBias);
        WIV_PB(z)=mean(WIV_Pbias);
        WIV_HB(z)=mean(WIV_Hbias);
        PLE_WIVP(z)=mean(PLE_WPbias);
        PLE_WIVH(z)=mean(PLE_WHbias);
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
Sigma = 0: 2: 180;
PB=[PLE_P;PBIAS;MPBIAS;PLE_WIVP;WIV_PB];%ple, bcple,ML. avplewiv, bcplvwiv
HB=[PLE_H;HBIAS;MHBIAS;PLE_WIVH;WIV_HB];
xlswrite('PB_rot', PB);        %位置误差保存
xlswrite('PH_rot', HB);        %角度误差保存
% 数据存储
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
plot(Sigma,PLE_P,'m*-',Sigma,PBIAS,'b-', Sigma, MPBIAS,'r.-', Sigma, WIV_PB,'g+--',Sigma, PLE_WIVP,'k--','linewidth',1.5)
%plot(Sigma,PLE_P,'m--',Sigma,PBIAS,'b-.', Sigma, MPBIAS,'r*--', Sigma, WIV_PB,'g+--','linewidth',1.5)
set(gca,'Fontsize',13);
h=legend('AVPLE','BCAVPLE','ML','BCAVPLE-WIV','AVPLE-WIV');
set(h,'Fontsize',12)
xlabel('Orientation of Unknown Node (degrees)');
ylabel('Location Error (m)');
axis([0 180 0 10]);
set(gca,'XTick',0:20:180);
grid on 

subplot(2,1,2)
plot(Sigma,PLE_H,'m*-',Sigma,HBIAS,'b-', Sigma, MHBIAS,'r.-', Sigma, WIV_HB,'g+--',Sigma, PLE_WIVH,'k--','linewidth',1.5)
%plot(Sigma,PLE_H,'m>--',Sigma,HBIAS,'b*--', Sigma, WIV_HB,'g+--',Sigma, PLE_WIVH,'k^--','linewidth',1.5)
set(gca,'Fontsize',13)
h=legend('AVPLE','BCAVPLE','ML','BCAVPLE-WIV','AVPLE-WIV');
set(h,'Fontsize',12)
xlabel('Orientation of Unknown Node (degrees)');
ylabel('Orienation Error (degrees)');
axis([0 180 0 10]);
set(gca,'XTick',0:20:180);
grid on 
