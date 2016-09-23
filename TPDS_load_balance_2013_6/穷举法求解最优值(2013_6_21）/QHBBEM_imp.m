clc,clear;                 %多元化电价下数据中心负载调度----穷举法
tic
format short               %线性规划问题
% %%%%%%转移延时%%%%%%%%%%%
% Pr1=[20.70 87.9 15.75;17.44 87.9 15.13;12.87 87.9 13.58;8.20 87.9 17.13;10.41 87.9 20.38;19.15 87.9	21.36;20.86 87.9 23.04;23.47 87.9 23.9;23.25 87.9 23.58;...
%     24.4 89 25;24.74 87.9 26.73;25.17 87.9 28.61;24.85 87.9 33.88;26.05 87.9 44.96;23.89 87.9 60.75;23.43 87.9 33.6;21.43 87.9 28.1;22.02 87.9 25.82;...
%     22.74 87.9 24.87;22.98 87.9 23.78;28.51 87.9 21.61;27.54 87.9 18.81;27.01 87.9 18.28;23.87 87.9 16.41];       % 表示电价input
Pr1=[34.8 33.2 28 25.6 25 24.6 22.8 24.8 27.2 31.85 33.4 35.2 36.2 36.2 36.8 37.2 38.5 40.4 41 42.4 41.7 47.1 45.2 36.6;...
74.73 74.73 74.73 74.73	74.73 74.73 74.73 74.73	74.73 74.73	74.73 74.73	74.73 74.73	74.73 74.73 74.73 74.73 74.73 74.73 74.73 74.73 74.73 74.73;...
19.35 18.2 17.87 17.99 20.54 21.95 25.79 29.22 34.38 36.58 38.23 37.93 39.5 49.94 55.45	50.93 42.49	36.91 36.73	36.09 30.4 27.14 25.73 21.03]';

%%%%%%转移延时%%%%%%%%%%%
DL1=[1.5 1 2 2.5 3 3.3 3.2 2 2.5 2.1 3.2 4.5 2.5 2.1 2 2 2.1 2.3 2.4 3 5.3 2 2 1.9;58 55 53 60 66 70 75 80 90 95 96 90 89 88 76 77 82 85 89 96 78 65 55 56;...
     35 32 30 40 50 55 60 72 78 82 87 89 82 75 64 65 62 67 88 75 55 48 44 35;52 56 50 46 54 61 68 70 76 85 90 95 88 84 77 76 74 76 78 86 72 65 60 57]/1000;
DL2=[58 55 53 60 66 70 75 80 90 95 96 90 89 88 76 77 82 85 89 96 78 65 55 56;2.5 1.8 2 2.3 2.5 2.8 2.7 3.2 3.4 3.6 3.8 4.4 3.5 3.2 2.4 2 2.2 2.5 2.9 2.8 3.6 1.5 1 1;...
     26 24 23 25 29 32 40 53 62 68 74 78 73 72 70 66 64 66 68 71 75 64 45 35;8 6 6.5 6 8 12 16 20 25 27 30 32 26 25 22 21 19 19 20 29 25 20 18 13]/1000;
DL3=[35 32 30 40 50 55 60 72 78 82 87 89 82 75 64 65 62 67 88 75 55 48 44 35;26 24 23 25 29 32 40 53 62 68 74 78 73 72 70 66 64 66 68 71 75 64 45 35;...
    1.1 1.2 1.5 2 2.4 2.8 2.5 3 3.6 3.7 2.9 2.8 2.7 2.6 2.1 2 2 2.5 3.4 1.5 1.6 1.5 1.1 1;38 38 36 37 35 46 56 65 73 75 78 75 72 74 70 64 67 69	78 73 64 55 46 40]/1000;
L1=30000*ones(1,24);
L2=40000*ones(1,24);
L3=40000*ones(1,24);
L4=30000*ones(1,24);
Fva=[];
for nj=1:24 %1,2,4
    u=20; 
    M=[5000 4000 4000];        %数据中心服务器最大数目  
    Pr=Pr1(nj,:);            %电价
    L=[L1(nj) L2(nj) L3(nj) L4(nj)];          %input
    n=length(M);                              %表示数据中心的个数
    cc=length(L);                              %表示前端服务器
    D=80/1000;                                    %延时上界 
    Tt=[DL1(:,nj)';DL2(:,nj)';DL3(:,nj)'];
    Tty=Tt;
    upr=unique(Pr);            %判断电价是否相同
    %DL=M*u;                    %服务器最大处理能力
    Po=133.99/10^6;                     %功率
    Aq1=[];
    Bq1=[];
    RU=inf;
    XH1=[];
    ki=0;
    %%%%%延时处理%%%%%%
    for i=1:n
        for j=1:cc
            if Tt(i,j)>=D
                ki=ki+1;
                Aq1(ki,:)=zeros(1,n+cc*n);
                Tt(i,j)=0;
                ij=(j-1)*n+i;
                XH1(ki)=ij;
                Aq1(ki,ij+n)=1;
                Bq1(ki)=0;
            end
        end
    end         %选择不能流向负载的数据线路
    T1=Tt(1,:);
    T2=Tt(2,:);
    T3=Tt(3,:);
    TT=Tt';
    Smt=Tt;
    YTT=Tt;
    f1=Pr*Po;               %目标函数多项式系数
    RU=f1*M';
    FMin=RU;
    FU=RU;
    %%%%%%%%%%%%%延时处理结束%%%%%%%%%%%
    for t1=1:cc
        for t2=1:cc
            for t3=1:cc
            StX=[];
            SUb=[];
            SLb=[];
            StF=[];
            Aeq=[];
            Beq=[];
            XR=[];   
            STm=[];
            A=[];
            B=[];
            A1=[];
            B1=[];
            A2=[];
            XH=[];
            Aqt=[];
            Bqt=[];
            AQ=[];
            BQ=[];
            XH2=[]; 
            M=[5000 4000 4000];        %数据中心服务器最大数目 
            L=[L1(nj) L2(nj) L3(nj) L4(nj)];          %input
            n=length(Pr);
            c=cc;
            Aqt=zeros(1,n+c*n);
            Bqt=0;
            jt=0;
            kt=0;
            for k1=1:c
               if T1(t1)<T1(k1)
                   jt=jt+1;
                   kt=kt+1;
                   kk=(k1-1)*n+1;  %k1表示前端服务器到数据中心1
                   XH2(kt)=kk;
                   Aqt(jt,kk+n)=1;%2n has been changed as n
                   Bqt(jt)=0;
               end
            end
            for k2=1:c
               if T2(t2)<T2(k2)
                   jt=jt+1;
                   kt=kt+1;
                   kk=(k2-1)*n+2;  %k2表示前端服务器
                   XH2(kt)=kk;
                   Aqt(jt,kk+n)=1;
                   Bqt(jt)=0;
               end
            end
            for k3=1:c
               if T3(t3)<T3(k3)
                   jt=jt+1;
                   kt=kt+1;
                   kk=(k3-1)*n+3;  %k3表示前端服务器到数据中心3
                   XH2(kt)=kk;
                   Aqt(jt,kk+n)=1;
                   Bqt(jt)=0;
               end
            end
            XH=[XH1,XH2];
            AQ=[Aq1;Aqt];
            BQ=[Bq1,Bqt];
            XH=sort(XH,2);
            Tm=[];
            Tm=[T1(t1),T2(t2),T3(t3)]; %选择的时间存储
            it=0;
            sg=0;
            ka=[];
            sm=[];
            H1=[];
            H2=[];
            H3=[];
            DM=[];
            for xj=1:c       %唯一流向的端口选择
                 y1=0;
                 gy=0;
                 SX=[];
                for xi=1:length(XH)
                   if XH(xi)<=n*xj && XH(xi)>n*(xj-1)
                      y1=y1+1;
                      SX(y1)=XH(xi);%负载不能流向的
                   end
                end
                if y1==(n-1)   %有唯一流向的
                    it=it+1;
                    for d=((xj-1)*n+1):n*xj
                        gy=0;
                        for y2=1:y1
                           if SX(y2)==d
                               gy=1;
                           end
                        end
                        if gy==0
                           dd=d;
                           dj=dd-(xj-1)*n;
                           H1(it)=dj;  %流向数据中心j
                           H2(it)=xj;  %端口
                           H3(it)=(xj-1)*n+1+n;
                       end
                    end
                    ka=[ka,xj];            %需要删除的端口i
                    Lm(it)=L(xj);          %负载存储
                    sg=1;
                end
            end
              LZ=zeros(n,1);
              for i=1:n
                  for j=1:length(H1)
                      if i==H1(j)
                         LZ(i)=LZ(i)+L(H2(j));  %同一端口负载存储
                      end
                  end     
              end  
              UH=unique(H1);
              for zi=1:length(UH)
                   zx=H1(zi);                                %数据中心
%                  dt=D-max(Smt(zx,:));
                   dt=D-Tm(zx); %%%%                        %最大的延时 为什么要最大？每次循环选择的时间中的最大dt=D-max(Smt(zx,:));
                   Ms=ceil((LZ(zx)+1/dt)/u);                 %打开的服务器数目
                   M(zx)=M(zx)-Ms;
                   DM(zi,:)=[zx,Ms];                        % 已经使用的数量存储
              end
              L([ka])=[];
             for hi=1:length(H3)
               hj=H3(hi)-(hi-1)*n;                           %删除端口i到数据中心的负载流
               AQ(:,[hj:hj+n-1])=[];
             end
             ai=0;
             Sq=[];
             if isempty(AQ)==1
                la=0;
             else
                la=length(AQ(:,1));
             end
             for i=1:la
                if ~all(AQ(i,:)==0)
                    ai=ai+1;
                    Sq(ai,:)=AQ(i,:);     %%error
                end
             end
             AQ=Sq;
             if ai>=1
               BQ=zeros(1,ai);
               else BQ=[];
             end   
    Aeq=[];
    Beq=[];  
    Tm=[];
    Tm=[T1(t1),T2(t2),T3(t3)];
    Tmm=Tm;
    if isempty(DM)==0
       Mb=unique(DM(:,1));
       for hn=1:length(Mb)  %已经有数据流
           Tm(hn)=-inf;
       end     
    end
    n=length(M);
    c=length(L);
    LL=[];
    for i=1:c
        Lp(1:n)=L(i);
        LL=[LL,Lp];
    end
    NM=find(M<0);
    if isempty(NM)==0
        c=-1;
    end
    if c>=1
    DL=M*u;
    f1=Pr*Po;               %目标函数多项式系数
    f2=zeros(1,c*n);
    f=[f1,f2];
    MOF=zeros(c,n);       
    Rmd1=[];
    xs=ones(1,n);           %产生n个1 
    for i1=1:c
        ZL=zeros(c,n);      %全是0阵列
        ZL(i1,:)=xs;    
        Rmd1=[Rmd1,ZL];
    end
    Aeq=[AQ;MOF,Rmd1];     %第一个约束条件: 负载守恒
    Beq=[BQ';L'];

    Dy=(D-Tm)';              %%%%%延时处理
    Rm2=eye(n,n);            %rmd 的系数
    Rmd2=[];
    for i=1:c
        Rmd2=[Rmd2,Rm2];
    end
    m2=-u*eye(n,n);          %m的系数
%     OF1=diag(1./Dy);         %ON/OFF的系数
%     A1=[m2,OF1,Rmd2];        %第二个约束条件：延时分析
      A=[m2,Rmd2];

%     OF2=-diag(DL);           %ON/OFF的系数
%     m3=zeros(n,n);
%     A2=[m3,OF2,Rmd2];        %第三个约束条件，附加约束

    %A=[A1;A2];
    %B=zeros(2*n,1);          %不等式右边
    D_bd=[];
    D_bd=1./Dy;               %delay constraint 不等式右边
    B=-D_bd;                 %不等式右边
    nl=c*n+n;
    Lb=[];
    Ub=[];
    Lb=zeros(nl,1);           %未知参数下界

    %OF3=ones(1,n);            %ON/OFF系数
    Ub=[M,LL]';
    nc=min(n,c);               %电价相同负载分配情况
    [x,fval,exitflag,ouput,lambda]=linprog(f,A,B,Aeq,Beq,Lb,Ub);
    else if c==0
        XR=zeros(n,1);
        c=cc;
        for i=1:n
            for j=1:length(DM(:,1)) 
                MD=DM(j,:);
                if i==MD(1)             %数据中心序号
                    XR(i)=XR(i)+MD(2);  %数据中心服务器数
                end
            end
        end
        %XR(n+1:n)=ones(n,1);% what does it mean?
        XR(n+1:n*c+n)=0;
        for ij=1:length(H1)
            i=H2(ij);                   %端口号
            j=H1(ij);                   %Data center
            nt=n+(i-1)*n+j;
            XR(nt)=Lm(ij);
        end
        XRS=XR;
        RU=Pr*XR(1:n)*Po;
        sg=0;
        exitflag=-1;
        else
            RU=inf;
            exitflag=-1;
        end
    end
    if  exitflag>0
        RL=fval;
        fp=fval;
        MRU=RU;
        k=1;
        ki=0;
        ni=0;
        StF(k)=MRU;
        [flag,S,H]=fNull(x,n);  
        hfg=1;
        SS=[];
        SH=[];
        SH=H;       %存储非整数的个数
        SS=S;       %存储节点的值
        UP=Ub;
        LP=Lb;
        if flag==0
            XRS=x;
            hfg=0;
        end
        if RL>=FMin;
            flag=0; hfg=0;
        end
     else flag=0; hfg=0;
    end
 while (flag|| hfg)
     UbP=[];
     LbP=[];
     UbP=Ub;                        %存储上一次
     LbP=Lb;
     BB=[];
     nh=length(H);  %%%%Srong Branch
     %fh=0;
%     for ih=1:nh
%         m=H(ih);                   % 确定第几个数
%         if m>=4 & m<=6
%            m=H(ih);
%            fh=1;
%            nh=0;
%            break;
%         end
%     end
    %if fh==0
       m=H(1);
       ih=1;
    %end
%      m=H(1);
%      ih=1;
     if hfg==1    %删除节点
       SH(ih)=[];
       SS(ih)=[];
     end
    for T=1:2                      %分枝结果
        if T==1
           Ub(m)=floor(S(ih));        %上界 
           BB(:,T)=Ub;
        else if T==2
          Lb(m)=floor(S(ih))+1;        %未知参数下界
          BB(:,T)=Lb;
        end
        end
    end
    
    for T=1:2
        if T==1
            Ub=BB(:,T);             %改变上界
            Lb=LbP;                 %下界不变
            [x,fval,exitflag,ouput,lambda]=linprog(f,A,B,Aeq,Beq,Lb,Ub);
            if exitflag>0
               x1=x;
               FF1=fval;
               [flag,S,H]=fNull(x1,n);
               n1=length(S);
               if (flag==0 && FF1<RU)
                   RU=floor(FF1*1000)/1000;
                   XR=x1;
               end
            else
               FF1=inf;
               n1=c*n+n;
               x1=[];
            end
        else if T==2
             Ub=UbP;
             Lb=BB(:,T);            %下界
             [x,fval,exitflag,ouput,lambda]=linprog(f,A,B,Aeq,Beq,Lb,Ub);
             if exitflag>0
               x2=x;
               FF2=fval;
               [flag,S,H]=fNull(x2,n);
               n2=length(S);
               if (flag==0 && FF2<RU)
                   RU=floor(FF2*1000)/1000;
                   XR=x2;
               end
             else
               FF2=inf;
               n2=c*n+n;
               x2=[];
             end
            end
        end
    end
    %% 获取最优整数
    if n1==0 && n2==0
        fz=0; flag=0;
     else fz=1;
    end
    F1=floor(FF1*1000)/1000;
    F2=floor(FF2*1000)/1000;
    Fmin=min(F1,F2);  %%%
    Fmax=max(F1,F2);
    %%%%%启发式函数
   if fz==1
   if (Fmin<RU && Fmin<FMin)
       fa=F1/(F1+F2);  %归一化
       fb=F2/(F1+F2);
       na=n1/(n1+n2);
       nb=n2/(n1+n2);
     if (fa+na)<(fb+nb)    %评价函数
        x=x1;
        Ub=BB(:,1);        %分支的选择
        Lb=LbP;
       if F2<RU
        StX(:,k)=x2;
        StF(k)=F2;         %原始节点的值
        SUb(:,k)=UbP;
        SLb(:,k)=BB(:,2);
        k=k+1;
       end
     else   
        x=x2;              %下界
        Lb=BB(:,2);
        Ub=UbP;
        if F1<RU
          StX(:,k)=x1;
          StF(k)=F1;
          SUb(:,k)=BB(:,1);
          SLb(:,k)=LbP;
          k=k+1;
        end
      end
    [flag,S,H]=fNull(x,n);
      else flag=0;
    end
   end
    if (flag==0)          %求出可行解
      [minF,ni]=min(StF);
       if (length(minF)>=1) && (RU>minF) && (minF<MRU)
          x=StX(:,ni);
          Ub=SUb(:,ni);
          Lb=SLb(:,ni);
          fp=minF;
          StF(ni)=[];
          StX(:,ni)=[];
          SUb(:,ni)=[];
          SLb(:,ni)=[];
          k=k-1;
       else x=[];
       end
    end
    [flag,S,H]=fNull(x,n);
    if (isempty(SH)~=1) && (isempty(XR)==1) && (isempty(x)==1)
        hfg=1;                         %判断节点搜索是否有效
        H=SH;
        S=SS;
        Ub=UP;
        Lb=LP;
      else hfg=0;
    end
    display('-----one------');
 end
   %%%%存储每次结果%%%%%ERROR
        display('----OK-----')
        FRU=0;
        if sg==1
          for ir=1:length(DM(:,1))
             Dm=DM(ir,:);   %服务器序号和数目大小
             FRU=Dm(2)*Po*Pr(Dm(1))+FRU;
          end
        end
        RU=RU+FRU;
        if RU<FMin && (isempty(XR)~=1)
               rmd=XR(n+1:n+c*n);
               SXR=XR;
               rx=[];
               rxt=[];
               for jt=0:(c-1)
                  for it=1:n
                    k2=jt*n+it;
                    rx(jt+1,it)=rmd(k2);
                  end
               end
%                 rxt=rx';
%                 FRU=0;
                if sg==1    %表明有负载单独流向数据中心
%                     for ir=1:length(DM(:,1))
%                         Dm=DM(ir,:);   %服务器序号和数目大小
%                         FRU=Dm(2)*Po*Pr(Dm(1))+FRU;
%                         XR(Dm(1))=Dm(2)+XR(Dm(1));
%                     end 
           %%%%%%%%%还原负载分配情况%%%%%%%%%%%
                    for kr=1:length(H1)   
                        jr=H2(kr);  %端口号
                        tr=H1(kr);  %Data center
                        Fz=zeros(1,n);
                        Fz(tr)=Lm(kr);
                        r=length(rx(:,1));
                        if jr<=r
                           rxj=rx([1:jr-1],:);
                           rxc=rx([jr:r],:);
                           rx=[rxj;Fz;rxc];
                           r=r+1;
                        else rx=[rx;Fz];
                        end
                    end
                  
                end
                rxt=rx';
             %% 求延时
                 fg=1;
                 Smt=YTT;
                 Tm=[];
                 for i=1:n
                     while (fg)
                          [tm,ri]=max(Smt');      %某一个数据中心最大延时的大小和位置
                          kk=ri(i);
                          if rxt(i,kk)>0.1  %最大值对应的rmd
                            Tm(i)=tm(i);       
                            fg=0;
                          else 
                              Smt(i,kk)=0;  %%%%
                               Ty=Smt(i,:);
                               if all(Ty<=0)
                                   Tm(i)=0;
                                   fg=0;
                               end
                          end     
                      end
                     fg=1;
                 end
                SMD=[];
                SMD=sum(rx);               %到某一个数据中心的负载数目
                for j=1:n
                    if SMD(j)>=1
                       Load(j)=SMD(j)+1./(D-Tm(j));
                    else Load(j)=0;
                    end
                end
               %%%%求打开的服务器数目%%%%%%%%
                Load=round(Load*100)/100;
                Server=ceil(Load/u);          %服务器数目
                RU=Server*Pr'*Po;
                if RU<FMin
                   Stm=Tmm;                    %存储时间 
                   FMin=RU;                
                   SSTm=Tm;                   %最优延时存储 
                   S_M=Server;                %服务器数目
                   RX=rxt;
                end

                
                    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     FRU=FRU+RU;
%                   if FRU<FMin
%                       UH=unique(H1);
%                      for i=1:length(UH);
%                          t_m=max(Smt(UH(i),:));
%                          dt1=(1/(D-t_m))/u; %tm开始的时候是最大；出现问题啦！
%                          dt2=(1/(D-Tmm(UH(i))))/u; %结束以后的tm
%                          dxr=ceil(dt2)-ceil(dt1);  %误差
%                          XR(UH(i))=XR(UH(i))+dxr;
%                      end
%                    FMin=FRU;                    
%                    SXR=XR;
%                    Stm=Tmm;
%                    RX=rx;
%                    Trx=rx';
%                   end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 end sg==1;
%                 if RU<FMin && (sg==0)
%                    FMin=RU;               
%                    Stm=Tmm;   %存储时间 
%                    SXR=XR;
%                    RX=rx;
%                    Trx=rx'; 
%                    display('-------Cong Get it-----');
%                 end

        end
              RU=FMin;    %非常重要的 important
              %%%清理空间%%%%
              XR=[];
              StX=[];
              SUb=[];
              SLb=[];
              StF=[]; 
          end
      end
   end

    %% 存储结果
    SPTM(nj,:)=SSTm;  %最优时间结果
    Tms(nj,:)=Stm;
    Fva=[Fva,FMin]; %函数值
    STM(nj,:)=S_M;  %服务器数目
 end
xlswrite('OPServer',STM);%%%server numbers
xlswrite('OPVal',Fva);  %optimal cost
xlswrite('OPTime',Tms);
xlswrite('SPtm',SPTM);  %optimal time
display('------*******the running is OK Now********---');
toc