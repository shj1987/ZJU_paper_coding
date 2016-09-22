%MSE and bias calculations for ML, OV-LS, OV-TLS and OV-CTLS estimators
%by Kutluyil Dogancay
%Feb 19, 2003
clear all
close all

%load simulation data
simdata;

r=zeros(N,2); %receiver position matrix
r_hat=zeros(N,2);
theta=zeros(N,1); %bearing vector
theta_hat=zeros(N,1);

rg=linspace(rxl,rxu,N);
j=0;
for i=rg
    j=j+1;
    ri=[i,a*i+b];
    r(j,:)=ri; %receiver position
    d=p'-ri; %range vector
    theta(j)=atan2(d(2),d(1)); %bearing angle
end

%plot of observer path and target location
figure(1);
plot(r(:,1),r(:,2),'^',p(1),p(2),'o');
h=gca; set(h,'FontSize',13);
axis(xy);
legend('Observer positions','Stationary target',2);
xlabel('x-axis (km)');
ylabel('y-axis (km)');

%no. of simulation runs
runs=10000;

%no. of simulations
sims=size(par,1);

%system matrix and data vector
A=zeros(N,2);
b=zeros(N,1);

%estimation results
pML=zeros(runs,2);
pLS=zeros(runs,2);
pTLS=zeros(runs,2);
pCTLS=zeros(runs,2);

ptrue=ones(runs,1)*p';

for i=1:sims
    i
    for j=1:runs
        %bearing errors (zero-mean white Gaussian)
        theta_hat=theta+par(i,1)*randn(N,1);
    
        %receiver position noise (zero-mean white Gaussian)
        r_hat=r+par(i,2)*randn(N,2);
        
        %system matrix with errors
        A=[sin(theta_hat),-cos(theta_hat)];

        %data vector with errors
        b=(r_hat(:,1).*sin(theta_hat))-(r_hat(:,2).*cos(theta_hat));

        %LS estimate of emitter location
        p_est=lscov(A,b);
        pLS(j,:)=p_est';
    
        %TLS estimate of emitter location
        p_est=tls(A,b);
        pTLS(j,:)=p_est';
    
        %CTLS 
        p_est=ctls(A,r_hat,mrec,errthr);
        pCTLS(j,:)=p_est';
    
        %ML
        %p_est=nlls(pinit,theta_hat,r_hat,par(i,1)^2*ones(N,1),par(i,3),par(i,4));
        p_est=fminsearch('mlcost',p,[],theta_hat,r_hat,ones(N,1));
        pML(j,:)=p_est';
    end

    [p_est,err]=ctls(A,r_hat,mrec,errthr);
    %[p_est,pp]=nlls(pinit,theta_hat,r_hat,par(i,1)^2*ones(N,1),par(i,3),par(i,4));
    %fminsearch('mlcost',pinit,[],theta_hat,r_hat,ones(N,1))
    
    figure;
    %subplot(2,1,1), plot(err);
    %subplot(2,1,2), plot(pp);
    plot(err); title('CTLS convergence');
    
    LSmean=mean(pLS);
    TLSmean=mean(pTLS);
    CTLSmean=mean(pCTLS);
    MLmean=mean(pML);

    LSbias(i)=norm(p-LSmean'); 
    TLSbias(i)=norm(p-TLSmean'); 
    CTLSbias(i)=norm(p-CTLSmean'); 
    MLbias(i)=norm(p-MLmean');
    
    LSMSE(i)=sum(sum((pLS-ptrue).^2,2))/runs;
    TLSMSE(i)=sum(sum((pTLS-ptrue).^2,2))/runs;
    CTLSMSE(i)=sum(sum((pCTLS-ptrue).^2,2))/runs;
    MLMSE(i)=sum(sum((pML-ptrue).^2,2))/runs;

end

Biases=[LSbias',TLSbias',CTLSbias',MLbias']

MSEs=[LSMSE',TLSMSE',CTLSMSE',MLMSE']

%save results to file 'simresults'
%save simresults Biases MSEs;

figure; 
plot(par(1:5,1)/pi*180,Biases(1:5,1),'x:',...
    par(1:5,1)/pi*180,Biases(1:5,2),'*-.',...
    par(1:5,1)/pi*180,Biases(1:5,3),'s--',...
    par(1:5,1)/pi*180,Biases(1:5,4),'o-')
h=gca; set(h,'FontSize',13);
xlabel('Bearing noise standard deviation (degrees)');
ylabel('Norm of estimation bias (km)');
legend('LS','TLS','CTLS','ML',2)

figure; 
plot(par(1:5,1)/pi*180,MSEs(1:5,1),'x:',...
    par(1:5,1)/pi*180,MSEs(1:5,2),'*-.',...
    par(1:5,1)/pi*180,MSEs(1:5,3),'s--',...
    par(1:5,1)/pi*180,MSEs(1:5,4),'o-')
h=gca; set(h,'FontSize',13);
xlabel('Bearing noise standard deviation (degrees)');
ylabel('MSE');
legend('LS','TLS','CTLS','ML',2)

figure; 
plot(par(6:11,2),Biases(6:11,1),'x:',...
    par(6:11,2),Biases(6:11,2),'*-.',...
    par(6:11,2),Biases(6:11,3),'s--',...
    par(6:11,2),Biases(6:11,4),'o-')
h=gca; set(h,'FontSize',13);
xlabel('Observer position error standard deviation (km)');
ylabel('Norm of estimation bias (km)');
legend('LS','TLS','CTLS','ML',0)

figure; 
plot(par(6:11,2),MSEs(6:11,1),'x:',...
    par(6:11,2),MSEs(6:11,2),'*-.',...
    par(6:11,2),MSEs(6:11,3),'s--',...
    par(6:11,2),MSEs(6:11,4),'o-')
h=gca; set(h,'FontSize',13);
xlabel('Observer position error standard deviation (km)');
ylabel('MSE');
legend('LS','TLS','CTLS','ML',2)

