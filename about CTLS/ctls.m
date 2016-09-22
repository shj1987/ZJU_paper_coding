function [x,err]=ctls(A,R,rec,thr);
%   X=CTLS(A,B) Constrained TLS estimate for emitter location problem
%   by Kutluyil Dogancay
%   December 19, 2002

%err=zeros(rec,1);
[M,N]=size(A);
rk=rank(A);
b=diag(A*R');

i=1;
flg=1;

while (i<=rec) & flg,
    C=[A,-b];
    [U,S,V]=svd(C);
    
    err(i)=S(rk+1,rk+1)/sqrt(sum(diag(S(1:rk,1:rk)).^2));
    if err(i)>thr,
        i=i+1;
        S(rk+1,rk+1)=0; 
        C0=U*S*V'; %reduced rank approximation (contains perturbed A and b)

        %project perturbed A onto the unit circle
        A=C0(:,1:2);
        L2=sqrt(sum(A.^2,2))*[1,1];
        A=A./L2;
        b=diag(A*R'); %projected b
    else flg=0; 
    end
end

%CTLS estimate
C=[A,-b];
[U,S,V]=svd(C);
x=V(1:N,rk+1)/V(N+1,rk+1); %CTLS estimate

%x=pinv(A)*b;

