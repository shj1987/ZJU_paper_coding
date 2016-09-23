clc;  
clear all
format long
close all
%%%%%信标数目对定位精度的影响
%% 设定全局变量
global SIGMA
%% 设定基本的参数，比如蒙特卡洛的个数，DOA误差范围，声源个数等参数
M =200;                         %Monte Caro times
%Unodes=[100,90,pi/3];            %未知节点的位置信息Unodes=[170,120,pi/3]; 
R=300;                            %感知半径
threshold = 150;
%% 在不同声源数目下进行定位
z=0;
SIGMA=2;
nm=0;
 for i=1:M    
        nm=nm+1;
        num=0;
        node=round(rand(1,2)*100); 
        hd=round(rand(1,1)*100*pi/2)/100;
        Unodes=[node hd];
        for BN = 10: 40 : 250 
            kk=1;
            ei=0;
            pi1=0;
            SPbias=[];
 %% %%%%%%%%%%% 未知节点对声源和信标节点的DOA测量值%%%%%%%%%%%%%%%%%%%%
            S_Bcon=round(rand(BN,2)*100); 
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
                [U,P_Bias, H_Bias] = PLS(SBStore,DOA,Unodes(1:2),Unodes(3));%使用Psedo-Linear Least Squares.
                [UBC,OV,BP_Bias, BH_Bias] = BCPLS(SBStore,DOA,Unodes(1:2),Unodes(3));%使用Bias Compensated PseudoLinear Least Squares.
                [IVP_Bias, IVH_Bias] = WIV(SBStore,DOA,Unodes(1:2),Unodes(3),UBC);
                [PLE_WIVP, PLE_WIVH] = PLEWIV(SBStore,DOA,Unodes(1:2),Unodes(3),U);%PLE-WIV
                [MP, MH]=BMLF(SBStore,DOA,Unodes(1:2),Unodes(3),OV);
                if  MP< threshold && MH< 10 && P_Bias< threshold
                    num= num + 1;
                    PLS_PBias(nm,num) = P_Bias;  %存储结果，运行次数num
                    PLS_HBias(nm,num) = H_Bias;
                    BCPL_PBias(nm,num)=BP_Bias;
                    BCPL_HBias(nm,num)=BH_Bias;
                    WIV_PBias(nm,num)=IVP_Bias;
                    WIV_HBias(nm,num)=IVH_Bias;
                    PLE_WP(nm,num)=PLE_WIVP;
                    PLE_WH(nm,num)=PLE_WIVH;
                    ML_PBias(nm,num)=MP;
                    ML_HBias(nm,num)=MH;
                else
                    nm=nm-1;
                    break;
                end
            end
            DOA=[];
            SBStore=[];
        end
        %display('------the programming is running now-----');
%         z=z+1;
%         PBIAS(z)= mean(PLS_PBias); %位置误差
%         HBIAS(z)= mean(PLS_HBias);
%         BC_PBIAS(z)= mean(BCPL_PBias);
%         BC_HBIAS(z)= mean(BCPL_HBias);
%         WIV_PBIAS(z)= mean(WIV_PBias);
%         WIV_HBIAS(z)= mean(WIV_HBias);
%         PLEWIV_PBIAS(z)= mean(PLE_WP);
%         PLEWIV_HBIAS(z)= mean(PLE_WH);
% %         ML_PBIAS(z)= mean(ML_PBias);
% %         ML_HBIAS(z)= mean(ML_HBias);
%         PLS_PBias=[];
%         PLS_HBias=[];
%         BCPL_PBias=[];
%         BCPL_HBias=[];
%         WIV_Pbias=[];
%         WIV_Hbias=[];
%         ML_HBias=[];
%         ML_PBias=[];
%         PLE_WP=[];
%         PLE_WH=[];
if mod(i,50)==0
    display('------the programming is running now-----');
end
 end
 PBIAS= mean(PLS_PBias); %位置误差
 HBIAS= mean(PLS_HBias);
 BC_PBIAS= mean(BCPL_PBias);
 BC_HBIAS= mean(BCPL_HBias);
 ML_PB=mean(ML_PBias)
 ML_HB=mean(ML_HBias);
 WIV_PBIAS= mean(WIV_PBias)
 WIV_HBIAS= mean(WIV_HBias);
 PLEWIV_PBIAS= mean(PLE_WP);
 PLEWIV_HBIAS= mean(PLE_WH);
BN =  10: 40 : 250;
%% 数据存储
PB=[PBIAS;BC_PBIAS;ML_PB;WIV_PBIAS; PLEWIV_PBIAS];%ple, bcple, ML, bcavplewiv, avplvwiv
HB=[HBIAS;BC_HBIAS;ML_HB;WIV_HBIAS; PLEWIV_HBIAS];
xlswrite('PBia',PB);  %PLE位置误差保存
xlswrite('HBia',HB);  %PLE角度误差保存

% xlswrite('PLE_PE',PBIAS);  %PLE位置误差保存
% xlswrite('PLE_HE',HBIAS);  %PLE角度误差保存
% xlswrite('BCPLE_PE',BC_PBIAS); %BCPLE位置误差保存
% xlswrite('BCPLE_HE',BC_HBIAS); %BCPLE角度误差保存
% xlswrite('BCPLEWIV_PE',WIV_PBIAS); %BCPLE位置误差保存
% xlswrite('BCPLEWIV_HE',WIV_HBIAS); %BCPLE角度误差保存
% xlswrite('PLEWIV_PE1',PLEWIV_PBIAS); %BCPLE位置误差保存
% xlswrite('PLEWIV_HE1',PLEWIV_HBIAS); %BCPLE角度误差保存
% xlswrite('ML_PE',ML_PBIAS); %BCPLE位置误差保存
% xlswrite('ML_HE',ML_HBIAS); %BCPLE角度误差保存


%% 图形显示2  节点方差
figure(1)
subplot(2,1,1)
plot(BN, PBIAS,'m>--',BN, BC_PBIAS,'b*--',BN, WIV_PBIAS,'go--',BN, PLEWIV_PBIAS,'kd--',BN, ML_PB,'rs--','linewidth',1.5)
%plot(BN, WIV_PBIAS,'g+--',BN, PLEWIV_PBIAS,'ro--','linewidth',1.5)
set(gca,'Fontsize',14);
legend('PLE','BCPLE','BCPLE-WIV','PLE-WIV','ML');
%legend('BCPLE-WIV','PLE-WIV');
xlabel('Number of Beacons（N)');
ylabel('Location Error (m)');
grid on 

subplot(2,1,2)
plot(BN, HBIAS,'m>--',BN, BC_HBIAS,'b*--',BN, WIV_HBIAS,'go--',BN, PLEWIV_HBIAS,'kd--',BN, ML_HB,'rs--','linewidth',1.5)
%plot(BN, WIV_HBIAS,'g+--',BN, PLEWIV_HBIAS,'ro--','linewidth',1.5)
set(gca,'Fontsize',14);
%legend('BCPLE-WIV','PLE-WIV');
legend('PLE','BCPLE','BCPLE-WIV','PLE-WIV','ML');
xlabel('Number of Beacons (N)');
ylabel('Orienation Error (degrees)');
grid on 
