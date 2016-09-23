clc;  
clear all
format long
close all
%%%%%信标数目对定位精度的影响
%% 设定全局变量
global SIGMA
%% 设定基本的参数，比如蒙特卡洛的个数，DOA误差范围，声源个数等参数
M =2000;                         %Monte Caro times
%Unodes=[100,90,pi/3];            %未知节点的位置信息Unodes=[170,120,pi/3]; 
R=300;                            %感知半径
threshold = 100;
%% 在不同声源数目下进行定位
z=0;
nm=0;
SIGMA=2;
BNN=410;
 for i=1:M
        num=0;                    
        nm=nm+1;
        node=round(rand(1,2)*100);
        Hdr=round(rand(1,1)*100*pi/2)/100;
        Unodes=[node,Hdr];
        for BN=10: 50 : BNN    
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
                [PLE_WIVP, PLE_WIVH] = WIV(SBStore,DOA,Unodes(1:2),Unodes(3),U);
                [MP, MH]=BMLF(SBStore,DOA,Unodes(1:2),Unodes(3),OV);
                if  P_Bias< threshold && MH< 10 && MP<threshold && PLE_WIVH<40 && IVH_Bias<40
                    num= num + 1;
                    %[TP_Bias,TH_Bias]=BMLThy_Node(SBStore,DOA,Unodes);
%                     CRLB_Pbias(nm,num)= TP_Bias;
%                     CRLB_Hbias(nm,num)= TH_Bias;
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
  if mod(i,100)==0
   display('------the programming is running now-----');
  end
end
% CRLB_P=mean(CRLB_Pbias);
% CRLB_H=mean(CRLB_Hbias);
PBIAS= mean(PLS_PBias); %位置误差
HBIAS= mean(PLS_HBias);
BC_PBIAS= mean(BCPL_PBias);
BC_HBIAS= mean(BCPL_HBias);
WIV_PBIAS= mean(WIV_PBias);
WIV_HBIAS= mean(WIV_HBias);
PLEWIV_PBIAS= mean(PLE_WP);
PLEWIV_HBIAS= mean(PLE_WH);
ML_PBIAS= mean(ML_PBias);
ML_HBIAS= mean(ML_HBias);
BN =  10: 50 : BNN;
%% 数据存储
PB=[PBIAS;BC_PBIAS;ML_PBIAS;WIV_PBIAS;PLEWIV_PBIAS];%ple, bcple, ML, avplewiv, bcplvwiv
HB=[HBIAS;BC_HBIAS;ML_HBIAS;WIV_HBIAS;PLEWIV_HBIAS];
xlswrite('PB_N', PB);        %位置误差保存
xlswrite('PH_N', HB);        %角度误差保存
% xlswrite('CRLB_PE',CRLB_P);  %PLE位置误差保存
% xlswrite('CRLB_HE',CRLB_H);  %PLE角度误差保存
% xlswrite('PLE_PE',PBIAS);  %PLE位置误差保存
% xlswrite('PLE_HE',HBIAS);  %PLE角度误差保存
% xlswrite('BCPLE_PE',BC_PBIAS); %BCPLE位置误差保存
% xlswrite('BCPLE_HE',BC_HBIAS); %BCPLE角度误差保存
% xlswrite('BCPLEWIV_PE',WIV_PBIAS); %BCPLE位置误差保存
% xlswrite('BCPLEWIV_HE',WIV_HBIAS); %BCPLE角度误差保存
% xlswrite('PLEWIV_PE',PLEWIV_PBIAS); %BCPLE位置误差保存
% xlswrite('PLEWIV_HE',PLEWIV_HBIAS); %BCPLE角度误差保存
% xlswrite('ML_PE',ML_PBIAS); %BCPLE位置误差保存
% xlswrite('ML_HE',ML_HBIAS); %BCPLE角度误差保存

%% 图形显示2  节点方差
figure(1)
subplot(2,1,1)
plot(BN, PBIAS,'m>--',BN, BC_PBIAS,'b*--',BN, WIV_PBIAS,'g+--',BN, PLEWIV_PBIAS,'kd--',BN, ML_PBIAS,'rs--','linewidth',1.5)
%plot(BN, PBIAS,'m>--',BN, BC_PBIAS,'b*--',BN, WIV_PBIAS,'go--',BN, PLEWIV_PBIAS,'kd--','linewidth',1.5)
set(gca,'Fontsize',13);
legend('PLE','BCPLE','BCPLE-WIV','PLE-WIV','ML');
%legend('PLE','BCPLE','BCPLE-WIV','PLE-WIV');
xlabel('The Number of Beacons (N)');
ylabel('Location Error (m)');
axis([10 410 0 4]);
set(gca,'XTick',10:50:410);
grid on 

subplot(2,1,2)
plot(BN, HBIAS,'m>--',BN, BC_HBIAS,'b*--',BN, WIV_HBIAS,'g+--',BN, PLEWIV_HBIAS,'kd--',BN, ML_HBIAS,'rs--','linewidth',1.5)
%plot(BN, HBIAS,'m>--',BN, BC_HBIAS,'b*--',BN, WIV_HBIAS,'go--',BN, PLEWIV_HBIAS,'kd--','linewidth',1.5)
set(gca,'Fontsize',13);
legend('PLE','BCPLE','BCPLE-WIV','PLE-WIV','ML');
xlabel('The Number of Beacons (N)');
ylabel('Orienation Error (degrees)');
axis([10 410 0 4]);
set(gca,'XTick',10:50:410);
grid on 
