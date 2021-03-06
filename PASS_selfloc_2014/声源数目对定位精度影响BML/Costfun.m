function J = Costfun(x,Lca,Agl)
% 最大似然估计的代价函数
% x：带估计的参数，1x3
% Lca：声源和信标节点的位置nx2
% Agl：测得的DOA值，1xn
    J = 0;
    C = [];
    for i = 1 : length(Agl)
        ARCTAN = atan2(Lca(i,2)-x(2),Lca(i,1)-x(1));
        C = Agl(i)-ARCTAN + x(3);  %弧度计算
        J = J + C^2;
    end