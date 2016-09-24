clc;  %%%利用静态声源辅助实现单个未知节点定位
clear all
format long
close all

%% 设定全局变量
global SIGMA
global M
%% 设定基本的参数，比如蒙特卡洛的个数，DOA误差范围，声源个数等参数
M = 1000;                        %Monte Caro times
NB1=[20 30 20 30];                           %信标节点
SS=[40 40 60 60];
R=100;                           %感知半径
threshold = 100;
err=3;
SIGMA=2;
%% 在不同声源数目下进行定位
 for ii=1:length(NB1)
    NB=NB1(ii);
    ASource=SS(ii);
    z=0;
    for Ss=20:10:80
        NU=Ss;
        num=0;
        nm=0;
        for i=1:M           
            SLocation=round(rand(ASource,2)*500); %多个声源位置分布500*500
            BL=round(rand(NB,2)*500);
            BA=round(rand(NB,1)*pi*10)/10;
            Beacon=[BL,BA];
            UL=round(rand(NU,2)*500);
            UA=round(rand(NU,1)*pi*10)/10;  
            Unode=[UL,UA];
            kk=1;
            ei=0;
            pi1=0;
            SPbias=[];
            S_Source=[];
            A_Source=[];
     %% 信标定位静态声源与未知节点%%%%%%%%
            for j=1: ASource
                for k=1:NB
                    BBeacon=Beacon(k,:);           %信标信息
                    SSource=SLocation(j,:);        %声源信息
                    if sqrt((SSource(1)-BBeacon(1))^2+(SSource(2)-BBeacon(2))^2)<=R   %感知半径范围之内
                       XITA= DOAMeasure(SSource,BBeacon); % 测声源的DOA值子函数
                       XXITA(kk)=XITA+BBeacon(3);             %绝对坐标系下的角度
                       Beacon_use(kk,:)= BBeacon;             %能够测到的信标节点
                       kk=kk+1;
                    end
                end
     %% 信标节点对静态声源进行定位
                 if (kk-1)>=2
                    OV1 = TLS_Single(XXITA,Beacon_use(:,1:2));    %声源估计位置存储
                    if all(OV1==0)~=1
                    Spt = BML1(XXITA,Beacon_use(:,1:2),OV1);
                    pc=sqrt((SSource(1)-Spt(1))^2+(SSource(2)-Spt(2))^2); %位置偏差
                    else pc=threshold+10;
                    end
                    if all(OV1==0)~=1 && pc<threshold  %设定阈值
                       ei=ei+1;
                       pi1=pi1+1;
                       SPbias(pi1)=pc;
                       S_Source(ei,:)=Spt;  %已经定位的声源存储
                       A_Source(ei,:)=SSource;%%准确位置
                    end
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
            %nm=length(AB_Source(:,1))
            cnt=0;
            for ij=1:NU
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
                if ku>=3                                
      %% %%%%%%%%%%%%%%%%%用已知的声源对未知节点进行定位BML   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    [P_Bias, H_Bias]= AVTLS(SBStore,DOA,Unodes(1:2),Unodes(3));
                    %[P_Bias, H_Bias] = BML(SBStore,DOA,Unodes(1:2),Unodes(3), OV);%使用Psedo-Linear Least Squares.
                    if  P_Bias< threshold && H_Bias<20
                        num= num + 1;
                        PLS_PBias(num) = P_Bias;  %存储结果，运行次数num
                        PLS_HBias(num) = H_Bias;
%                         PLE_P(num)=PB;
%                         PLE_H(num)=HB;
                        Snum(num)=ku;
                        cnt=cnt+1;
                    end
                end  
                DOA=[];
                SBStore=[];
                S_Source=[];  %一次结束
            end
        end
        display('---------the programing is running now-----------------');
        S_Error=[];
        z=z+1;
        PBIAS(ii,z) = mean(PLS_PBias); %位置误差
        HBIAS(ii,z)= mean(PLS_HBias);
        SBIAs(ii,z)=mean(SPavg);       %声源平均偏差
        SPavg=[];
        PLS_PBias=[];
        PLS_HBias=[];
        PLE_P=[];
        PLE_H=[];
        Snum=[];
        PI=[];
    end   
 end
Sigma=20:10:80;
%% 数据存储
xlswrite('Beacon_PE',PBIAS);%位置误差保存
xlswrite('Beacon_HE',HBIAS);%角度误差保存

%% 图形显示1  位置节点平均误差
%% 图形显示2  节点方差
figure(1)
subplot(2,1,1)
plot(Sigma,PBIAS(1,:),'b*--', Sigma,PBIAS(2,:),'kd--', Sigma,PBIAS(3,:),'r^--',Sigma,PBIAS(4,:),'go--','linewidth',1.5)
%plot(Sigma,PBIAS(1,:),'b*--', Sigma,PBIAS(2,:),'rd--',Sigma,PBIAS(3,:),'go--','linewidth',1.5)
set(gca,'Fontsize',13)
%legend('K=9','K=13','K=15','K=17','K=18');
legend('N=20,K=40','N=30,K=40','N=20,K=60','N=30,K=60');
xlabel('The Number of Unknown Nodes');
ylabel('Location Error (m)');
grid on 

subplot(2,1,2)
plot(Sigma, HBIAS(1,:),'b*--', Sigma, HBIAS(2,:),'kd--', Sigma, HBIAS(3,:),'r^--',Sigma, HBIAS(4,:),'go--','linewidth',1.5)
%plot(Sigma, HBIAS(1,:),'b*--', Sigma, HBIAS(2,:),'rd--',Sigma, HBIAS(3,:),'go--','linewidth',1.5)
set(gca,'Fontsize',13)
%legend('K=9','K=13','K=15','K=17','K=18');
legend('N=20,K=40','N=30,K=40','N=20,K=60','N=30,K=60');
xlabel('The Number of Unknown Nodes');
ylabel('Orienation Error (degrees)')
grid on 



