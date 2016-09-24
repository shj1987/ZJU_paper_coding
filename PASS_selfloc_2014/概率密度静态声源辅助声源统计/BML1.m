function  BMLEst = BML1(DOAM,Source,BMLInit)
f = @(x)Costfun1(x, Source,DOAM); %´ú¼Ûº¯Êý
options = optimset('Display' , 'off' , 'TolFun' , 1e-6);
BMLEst = fminsearch(f, BMLInit, options);

