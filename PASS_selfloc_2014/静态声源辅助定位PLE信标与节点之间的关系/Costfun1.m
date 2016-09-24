function J = Costfun1(x,Lca,ADOA)
% 最大似然估计的代价函数
% x：带估计的参数，1x3
% Lca：声源和信标节点的位置nx2
% ADOA：测得的DOA值，1xn
  J = 0;
  for i = 1 : length(ADOA)
      ARCTAN = atan2(x(2)-Lca(i,2),x(1)-Lca(i,1));
      C = ADOA(i)-ARCTAN;  %弧度计算
      J = J + C^2;
  end