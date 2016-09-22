function [x,P,sv]=tls(A,b);
%   [x,P,sv]=TLS(A,b) Total least squares estimate of x in (A+D)x=b+d
%   Kutluyil Dogancay 
%   December 18, 2002
%   Nov 12, 2005 (svd ratio)
[M,N]=size(A);
B=[A,b];
if rank(B)<=N, disp('rank deficient matrix'); end
%R=rank(A);
[U,S,V]=svd(B,0);
x=-V(1:N,N+1)/V(N+1,N+1);    %TLS solution
P=-S(N+1,N+1)*U(:,N+1)*V(:,N+1)';   %TLS Perturbation matrix
sv=S(N+1,N+1);  %min-to-max SVD ratio for [A,b]

