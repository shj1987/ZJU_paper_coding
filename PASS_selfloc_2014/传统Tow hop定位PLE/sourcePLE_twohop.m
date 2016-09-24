clc;  %%%利用静态声源辅助实现单个未知节点定位
clear all
format long
close all

%% 设定全局变量
display('----start now-----');
global SIGMA
global M
%% 设定基本的参数，比如蒙特卡洛的个数，DOA误差范围，声源个数等参数
M = 1000;                %Monte Caro times
NumSource =50;            %声源的个数
NumB=30;                  %信标节点
R=100;                    %感知半径
threshold = 150;
err=4;
%% 在不同声源数目下进行定位
 for ii=1:length(NumSource)
    ASource=NumSource(ii);
    z=0;
    for SIGMA = 3: 1 : 3
        num=0;
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
            jj=1;
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
                       %DOA(kk)=XITA+BBeacon(3); 
                       Beacon_use(kk,:)= BBeacon;          %能够测到的信标节点
                       kk=kk+1;
                    end
                end
     %% 信标节点对静态声源进行定位
                 if (kk-1)>=3
                    %Spt = TLS_Single(DOA,Beacon_use(:,1:2));    %声源估计位置存储
                    [SPT,pc]=PLS1(Beacon_use(:,1:2),XXITA,SSource(1:2),SSource(3));
                    %pcer=sqrt((SSource(1)-Spt(1))^2+(SSource(2)-Spt(2))^2) %位置偏差
                    if pc(1)<threshold  %设定阈值
                       ei=ei+1;
                       pi1=pi1+1;
                       SPbias(pi1)=pc(1);
                       S_Source(ei,:)=SPT;  %已经定位的声源存储
                       A_Source(ei,:)=SSource(1:2);%%准确位置
                       cnt=cnt+1;
                       num= num + 1;
                       PLS_PBias(num) = pc(1);  %存储结果，运行次数num
                       PLS_HBias(num) = pc(2);
                    else Unode(jj,:)=SSource; jj=jj+1;
                    end
                 else
                     Unode(jj,:)=SSource;
                     jj=jj+1;
                 end
                 XXITA=[];
                 Beacon_use=[];   
                 kk=1;             
            end
            if isempty(SPbias)~=1
                nm=nm+1;
                SPavg(nm)=mean(SPbias);  %一次操作声源平均偏差            
            end
 %% %%%%%%%%%%% 未知节点对声源和信标节点的DOA测量值%%%%%%%%%%%%%%%%%%%%
            S_Bcon=Beacon(:,1:2);
            B_Source=[];
            AB_Source=[];
            B_Source=[S_Source;S_Bcon];         %信标和声源位置信息
            %nn=length(B_Source(:,1))
            AB_Source=[A_Source;S_Bcon];      %准确位置  
            if isempty(Unode)==1
                NumU=0;
            else
               NumU=length(Unode(:,1));
            end
%             for ij=1:NumU
%                 ku=0;
%                 for iu=1:length(B_Source(:,1))
%                      Unodes=[];
%                      BSL=AB_Source(iu,:);           %已定位的声源和信标
%                      E_BSL=B_Source(iu,:);          %估计位置
%                      Unodes=Unode(ij,:);
%                      if sqrt((BSL(1)-Unodes(1))^2+(BSL(2)-Unodes(2))^2)<=R                       
%                            ku=ku+1;
%                            DOA(ku)= DOAMeasure(BSL,Unodes);   % 测得的带有噪声DOA值
%                            SBStore(ku,:)=E_BSL;    %声源估计位置存储
%                      end
%                 end
%                 if ku>=3                                
%       %% %%%%%%%%%%%%%%%%%用已知的声源对未知节点进行定位BML   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     [P_Bias, H_Bias] = PLS(SBStore,DOA,Unodes(1:2),Unodes(3));%使用Psedo-Linear Least Squares.
%                     if  P_Bias< threshold 
%                         num= num + 1;
%                         PLS_PBias(num) = P_Bias;  %存储结果，运行次数num
%                         PLS_HBias(num) = H_Bias;
%                         Snum(num)=ku;
%                         cnt=cnt+1;
%                     end
%                 end  
%                 DOA=[];
%                 SBStore=[];
%                 S_Source=[];  %一次结束
%             end
            PI(i)=cnt/NumSource*100;
        end
        display('---------the programing is running now-----------------');
        S_Error=[];
        z=z+1;
%         Num(ii,z)=mean(Snum); %参与定位个数
        Prt(ii,z)=mean(PI);  %未知节点定位百分比覆盖率
        xlswrite('Port',Prt);
        PBIAS(ii,z) = mean(PLS_PBias); %位置误差
        xlswrite('PBIAS',PBIAS);
        PVAR(ii,z) =std(PLS_PBias);    %位置方差
        HBIAS(ii,z)= mean(PLS_HBias);
        HVAR(ii,z) =std(PLS_HBias);    %角度方差
        xlswrite('HBIAS',HBIAS);
        SBIAs(ii,z)=mean(SPavg);       %声源平均偏差
        SPavg=[];
        PLS_PBias=[];
        PLS_HBias=[];
        Snum=[];
        PI=[];
    end
    
 end
Sigma=1:err;
%% 数据存储
xlswrite('TwohopPLE_PE32',PBIAS);%位置误差保存
xlswrite('TwohopPLE_HE32',HBIAS);%角度误差保存
xlswrite('TwohopPLE_PRT32',Prt);

%% 图形显示1  位置节点平均误差
% figure(1)
% subplot(2,1,1)
% %plot(Sigma,PBIAS(1,:),'b*--', Sigma,PBIAS(2,:),'kd--', Sigma,PBIAS(3,:),'r^--',Sigma,PBIAS(4,:),'go--','linewidth',1.5)
% plot(Sigma,PBIAS(1,:),'b*--','linewidth',1.5)
% set(gca,'Fontsize',14)
% %legend('K=9','K=13','K=15','K=17','K=18');
% %legend('K=30','K=40','K=50');
% xlabel('Bearing Noise Standard Deviation (degree)');
% ylabel('Location Error (m)');
% xlim([1 err])
% grid on 
% % 
% subplot(2,1,2)
% %plot(Sigma, HBIAS(1,:),'b*--', Sigma, HBIAS(2,:),'kd--', Sigma, HBIAS(3,:),'r^--',Sigma, HBIAS(4,:),'go--','linewidth',1.5)
% plot(Sigma, HBIAS(1,:),'b*--','linewidth',1.5)
% set(gca,'Fontsize',14)
% %legend('K=9','K=13','K=15','K=17','K=18');
% %legend('K=30','K=40','K=50');
% xlabel('Bearing Noise Standard Deviation (degree)');
% ylabel('Orienation Angle Error (degree)')
% xlim([1 err])
% grid on 
display('-----OK***********NOW----')
% %%
% figure(2)
% plot(Sigma,Num(1,:),'b*--','linewidth',1.5)
% set(gca,'Fontsize',14)
% %legend('K=10','K=20','K=50');
% xlabel('Bearing Noise Standard Deviation (degree)');
% ylabel('The Number of Beacons and Nodes ');
% xlim([1 err])
% grid on 
% 
% %% 未知节点覆盖比例
% figure(3)
% plot(Sigma,Prt(1,:),'b*--','linewidth',1.5)
% set(gca,'Fontsize',14)
% %legend('K=10','K=20','K=50');
% xlabel('Bearing Noise Standard Deviation (degree)');
% ylabel('The Proportion of Effective Self_localization (%)');
% xlim([1 err])
% grid on 
% %%  声源的平均误差
% figure(3)
% plot(Sigma, SBIAs(1,:),'b*--','linewidth',1.6)
% set(gca,'Fontsize',17)
% xlabel('DOA误差，单位：°')
% ylabel('位置误差，单位：m')
% xlim([1 err])
% grid on 
display('-----OK***********NOW----')


