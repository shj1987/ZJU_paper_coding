%% 生成n个阵列观测到的m个声源的角度
function XITA = DOAMeasure(Target, Array)  %source 
    global SIGMA
    xita = atan2(Target(2)-Array(2), Target(1)-Array(1)); %绝对坐标系的值在[-180，180]之间根据象限确定
    Sigma = normrnd(0, SIGMA)*pi/180;  %产生标准差为SIGMA的角度误差
    XITA = xita - Array(3)+ Sigma;
  
