function [flag,S,H]=fNull(x,n)
j=1;
S=[];
H=[];
for i=(1+n):length(x)
    if abs(x(i)-round(x(i)))>0.001  %这个效果不错0.00001
        S(j)=x(i);           %非整数保存
        H(j)=i;              %标号
        j=j+1;
    end
end
if isempty(S)==1
    flag=0;
else 
    flag=1;
end