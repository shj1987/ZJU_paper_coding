clc,clear
global SIGMA
M = 1000;                        %Monte Caro times
NumSource =6;                    %声源的个数
%Beacon=[10,10,pi/6;170,40,pi/5;90,110,pi/3;15,180,pi/4]; %Beacon=[10,10,pi/6;170,40,pi/5;90,110,pi/3;15,180,pi/4];
Beacon1=round(rand(100,2)*200);
Beacon2=round(rand(100,1)*pi*10)/10;
Beacon=[Beacon1 Beacon2];
Bl=length(Beacon(:,1));         %信标个数 Beacon=[10,20,pi/3;0,160,pi/6;170,190,pi/5;130,10,pi/4;60,100,pi/5];
Unodes=[170,120,pi/3];            %未知节点的位置信息Unodes=[170,120,pi/3]; 
R=100;                            %感知半径
threshold = 150;
err=6;
%% 在不同声源数目下进行定位
for ii=1:length(NumSource)
    ASource=NumSource(ii);
    z=0;
    for SIGMA = 1: 1 : err
        num=0;
        nm=0;
        for i=1:M    
            SLocation=round(rand(ASource,2)*200); %多个声源位置分布200*200
            kk=1;
            ei=0;
            pi1=0;
            SPbias=[];
            Th_Err=[];
            S_Source=[];
            A_Source=[];
     %% 信标定位静态声源与未知节点%%%%%%%%
            for j=1: ASource
                for k=1:Bl
                    BBeacon=Beacon(k,:);      %信标信息
                    SSource=SLocation(j,:);        %声源信息
                    if sqrt((SSource(1)-BBeacon(1))^2+(SSource(2)-BBeacon(2))^2)<=R   %感知半径范围之内
                       [XITA,Accu_Ang]= Current_Radian(SSource,BBeacon); % 测声源的DOA值子函数
                       XXITA(kk)=XITA+BBeacon(3);         %绝对坐标系下的角度
                       AXXITA(kk)=Accu_Ang;                   %准确的角度
                       Beacon_use(kk,:)= BBeacon;             %能够测到的信标节点位置
                       kk=kk+1;
                    end
                end
     %% 信标节点对静态声源进行定位
                if (kk-1)>=2
                    Spt = TLS_Single(XXITA,Beacon_use(:,1:2));            %声源位置估计子函数   %声源估计位置存储
                    pc=sqrt((SSource(1)-Spt(1))^2+(SSource(2)-Spt(2))^2); %位置偏差
                    if all(Spt==0)~=1 && pc<threshold  %设定阈值
                       ei=ei+1;
                       pi1=pi1+1;
                       SPbias(pi1)=pc;
                       TErr= Expect_Err(AXXITA,Beacon_use(:,1:2),SSource);   %理论误差
                       Th_Err(pi1)=TErr;
                       S_Source(ei,:)=Spt;  %已经定位的声源存储
                       A_Source(ei,:)=SSource;%%准确位置
                    end
                end
                XXITA=[];
                AXXITA=[];
                Beacon_use=[];   
                kk=1;             
            end
            if isempty(SPbias)~=1
                nm=nm+1;
                SPavg(nm)=mean(SPbias);  %一次操作声源平均偏差   
                STErr(nm)=mean(TErr);
            end
               
        end
        display('----the programing is running now------');
        z=z+1;
        Real_Err(z)=mean(SPavg);
        Theo_Err(z)=mean(STErr);   
    end 
end
 Sig=1:err;
 xlswrite('Real_E',Real_Err);
 xlswrite('Theo_E',Theo_Err);
 figure(1)
 plot(Sig,Real_Err,'b*--', Sig,Theo_Err,'r^--','linewidth',1.5)
 set(gca,'Fontsize',14)
 legend('Real Error','Theory Error');
 xlabel('Bearing Noise Standard Deviation (degree)');
 ylabel('Location Error (m)');
 xlim([1 err])
 grid on 
 display('----OK Now------');
 
 