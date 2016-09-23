Beacon=[45,55;30,70;95,50;100,50;85,35;40,40;50,15;60,40;20,30;75,60];
Unodes=[100,90,pi/3];
plot(Beacon(:,1),Beacon(:,2),'bo','MarkerEdgeColor','b','MarkerFaceColor','g','MarkerSize',9);
hold on
plot(Unodes(1),Unodes(2),'r>','MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',9);
grid on
legend('Beacon','Unknown Node');
xlim([10 130]);
ylim([10 130]);
set(gca,'Fontsize',13);
xlabel('Horizonal Location of the Nodes (m)');
ylabel('Vertical Location of the Nodes (m)');
