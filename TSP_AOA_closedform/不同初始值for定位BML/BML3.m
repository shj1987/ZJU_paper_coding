function [BMLPbias, BMLHbias] = BML3(Source,DOAM, Array, Heading)
global SIGMA
if SIGMA<=4
BMLInit = [55,95,pi/5];
else if SIGMA==5
       BMLInit = [5,90,pi/5];
    else if SIGMA==6
             BMLInit = [20,15,pi/5]; 
        else if SIGMA==7
                BMLInit=[5,5,pi/15]; 
            else if SIGMA==8
                     BMLInit=[10,5,pi/15];
                else if SIGMA==9
                        BMLInit=[25,0,pi/20];
                   else if SIGMA==10
                            BMLInit=[45,10,pi/15];
                       end
                    end
                end
            end
        end
    end
end
f = @(x)Costfun(x,Source,DOAM); %´ú¼Ûº¯Êý
options = optimset('Display' , 'off' , 'TolFun' , 1e-5);
BMLRE = fminsearch(f, BMLInit, options);
BMLPbias = sqrt((BMLRE(1)-Array(1))^2+(BMLRE(2)-Array(2))^2);
H_Bias = (BMLRE(3)-Heading)*180/pi;
BMLHbias=abs(H_Bias);