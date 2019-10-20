function [  ] = test_filters_noBias(  )

addpath('..\rotation3d','..\Filters Without Bias','..\Generate Motions');

N = 60;
sf = 200;

parfor n = 1:N
    [qTrue,gyro] = genTrig(60,sf);
    qMea = genMea(qTrue);
    q0 = mulQua(qTrue(1,:),expQua([pi,0,0]));
    qEstM = MEKF(gyro,sf,q0,qMea);
    qEstE = KFEul(gyro,sf,q0,qMea);
    qEstQ = KFQua(gyro,sf,q0,qMea);
    qEstC = Comp(gyro,sf,q0,qMea);
    
    parsave(n,qTrue,gyro,qMea,qEstM,qEstE,qEstQ,qEstC,true);
end

parfor n = 1:N
    [qTrue,gyro] = genRotAxis(60,sf);
    qMea = genMea(qTrue);
    q0 = mulQua(qTrue(1,:),expQua([pi,0,0]));
    qEstM = MEKF(gyro,sf,q0,qMea);
    qEstE = KFEul(gyro,sf,q0,qMea);
    qEstQ = KFQua(gyro,sf,q0,qMea);
    qEstC = Comp(gyro,sf,q0,qMea);
    
    parsave(n,qTrue,gyro,qMea,qEstM,qEstE,qEstQ,qEstC,false);
end

rmpath('..\rotation3d','..\Filters Without Bias','..\Generate Motions');

end


function parsave(n,qTrue,gyro,qMea,qEstM,qEstE,qEstQ,qEstC,t)

if t
    save(strcat('C:\result-filterComp\filters_noBias\trig\',num2str(n),'.mat'),...
        'qTrue','gyro','qMea','qEstM','qEstE','qEstQ','qEstC');
else
    save(strcat('C:\result-filterComp\filters_noBias\rotAxis\',num2str(n),'.mat'),...
        'qTrue','gyro','qMea','qEstM','qEstE','qEstQ','qEstC');
end

end

