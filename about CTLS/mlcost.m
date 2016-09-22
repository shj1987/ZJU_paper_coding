function j=mlcost(p,b,r,var);
%   J=MLCOST(P,B,R,VAR) Maximum likelihood cost function
%   P is 2x1 location vector
%   B is the Nx1 bearing measurements vector
%   R is the Nx2 observer position matrix
%   VAR is the Nx1 bearing error variance vector
j=sum(((b-atan2(p(2)-r(:,2),p(1)-r(:,1))).^2)./var);
