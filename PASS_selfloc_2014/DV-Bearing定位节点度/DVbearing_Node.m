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
NumSource =[60 60 60 60 60];            %声源的个数
NumB=30;                  %信标节点
R=100;                    %感知半径
threshold = 100;
err=3;
%% 在不同声源数目下进行定位
 for ii=1:length(NumSource)
    ASource=NumSource(ii);
    z=0;
    for SIGMA = 2: 1: 2
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
                RBeacon=Beacon;
                for k=1:NumB
                    BBeacon=Beacon(k,:);           %信标信息
                    SSource=SLocation(j,:);        %node 信息
                    if sqrt((SSource(1)-BBeacon(1))^2+(SSource(2)-BBeacon(2))^2)<=R   %感知半径范围之内
                       XITA= DOAMeasure(BBeacon,SSource);   % 测声源的DOA值子函数
                       XXITA(kk)=XITA;                      %DOA
                       Beacon_use(kk,:)= BBeacon;           %能够测到的信标节点
                       dk=k-kk+1;
                       RBeacon(dk,:)=[];
                       kk=kk+1;
                    end
                end
                ReSource=SLocation;
                ReSource(j,:)=[];              %剩余的节点
                SSource=SLocation(j,:);        %node 信息
                nn=0;
                NNode=[];
                for k1=1:(NumSource-1)
                    RNode=ReSource(k1,:);
                    if sqrt((SSource(1)-RNode(1))^2+(SSource(2)-RNode(2))^2)<=R   %感知半径范围之内
                        nn=nn+1;
                        NNode(nn,:)=RNode;   %other nodes
                    end
                end
                if isempty(NNode)==1
                    Nl=0;
                else
                    Nl=length(NNode(:,1));
                end
                for k2=1:length(RBeacon(:,1))
                    nb=0;
                    for k3=1:Nl
                        RBn=RBeacon(k2,:);   %Beacon
                        Nod=NNode(k3,:);     %Node
                        if sqrt((RBn(1)-Nod(1))^2+(RBn(2)-Nod(2))^2)<=R   %感知半径范围之内
                           nb=nb+1;
                           SNode(nb,:)=Nod; %节点存储
                        end
                    end
                    if nb>=2
                       Thta=DOA2_Measure(RBn,SNode(1:2,:),SSource);
                       if isempty(Thta)~=1
                          Beacon_use(kk,:)= RBn; 
                          XXITA(kk)=Thta;
                          kk=kk+1;
                       end
                    end
                    SNode=[];
                end                     
                 ei=ei+1;
                 De(ei)=kk-1;
                 XXITA=[];
                 Beacon_use=[];   
                 kk=1;             
            end
            Deg(i)=mean(De);
            De=[];
        end
        display('---------the programing is running now-----------------');
        z=z+1;
        Degree(ii,z)=mean(Deg);  %未知节点定位百分比覆盖率
        PLS_PBias=[];
        PLS_HBias=[];
        Snum=[];
        Deg=[];
    end
    
 end
%% 数据存储
% xlswrite('DVBearing_PE',PBIAS);%位置误差保存
% xlswrite('DVBearing_HE',HBIAS);%角度误差保存
xlswrite('DV_Degree_30',Degree);

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
Degree
display('-----OK***********NOW----');


