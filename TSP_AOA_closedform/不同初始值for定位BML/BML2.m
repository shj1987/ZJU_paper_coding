function [BMLPbias, BMLHbias] = BML2(Source,DOAM, Array, Heading)
global SIGMA
if SIGMA<=4
BMLInit = [25,95,pi/5];
else if SIGMA==5
       BMLInit =  [5,10,pi/15];
    else if SIGMA==6
             BMLInit = [5,5,pi/15];
        else if SIGMA==7
                BMLInit=[5,5,pi/15];
            else if SIGMA==8
                     BMLInit=[10,5,pi/15]; 
                else if SIGMA==9
                       BMLInit=[5,0,pi/10];
                   else if SIGMA==10
                            BMLInit=[0,0,pi/10];
                       end
                    end
                end
            end
        end
    end
end

                                
f = @(x)Costfun(x,Source,DOAM); %´ú¼Ûº¯Êý
options = optimset('Display' , 'off' , 'TolFun' , 1e-5);
BMLEst = fminsearch(f, BMLInit, options);
BMLPbias = sqrt((BMLEst(1)-Array(1))^2+(BMLEst(2)-Array(2))^2);
H_Bias = (BMLEst(3)-Heading)*180/pi;
BMLHbias=abs(H_Bias);