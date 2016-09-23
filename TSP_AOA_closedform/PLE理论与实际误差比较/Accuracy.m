clc,clear;
PE=xlsread('PLE_PEall');
%Num=[5 7 9 11 13 15 16 17 18 19 20 21 22 23 24 25];
Num=5:2:25;
plot(Num, PE(:,1),'b*--', Num, PE(:,2),'kd--', Num, PE(:,3),'r^--',Num, PE(:,4),'go--',Num, PE(:,5),'mp--',Num, PE(:,6),'ys--','linewidth',1.5)
set(gca,'Fontsize',14)
legend('Bearing Noise=1бу','K=11','K=13','K=15');
xlabel('Bearing Noise Standard Deviation (degree)');
ylabel('Orienation Angle Error (degree)')
grid on 