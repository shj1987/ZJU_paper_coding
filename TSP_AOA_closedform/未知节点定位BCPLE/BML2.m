function [BMLPbias, BMLHbias] = BML2(Source,DOAM, Array, Heading)
BMLInit2 = [40,25,pi/5];
f = @(x)Costfun(x,Source,DOAM); %´ú¼Ûº¯Êý
options = optimset('Display' , 'off' , 'TolFun' , 1e-6);
BMLRE = fminsearch(f, BMLInit2, options);
BMLPbias = sqrt((BMLRE(1)-Array(1))^2+(BMLRE(2)-Array(2))^2);
H_Bias = (BMLRE(3)-Heading)*180/pi;
BMLHbias=abs(H_Bias);