function [flag,S,H]=fNull(x,n)
j=1;
S=[];
H=[];
%整数精确位数；根据数据中心数目可以调整，在0.0001-0.01之间；
%当数据中心数目n较少时，参考为0.001-0.0001;当数据中心数目增加时，参考为0.01-0.001;
factor=0.01;  
for i=(1+n):length(x)
    if abs(x(i)-round(x(i)))>factor %0.001可以的
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