function [  ] = test_vecMea_noBias(  )

addpath('..\rotation3d','..\Filters Without Bias',...
    '..\Generate Motions','..\Vector to Rotation');

% global settings
N = 60;
sf = 200;
tTot = 60;
vRef = [0,1;0,0;1,0];

parfor n = 1:N
    % true trajectory and measurement
    [qTrue,gyro] = genTrig(tTot,sf);
    vMea1 = genMea(qTrue,'v',vRef(:,1),'vstd',0.1);
    vMea2 = genMea(qTrue,'v',vRef(:,2),'vstd',0.1);
    vMea = cat(2,reshape(vMea1',3,1,[]),reshape(vMea2',3,1,[]));
    
    % filters
    q0 = mulQua(qTrue(1,:),expQua([pi,0,0]));
    qEstM = MEKF(gyro,sf,q0,[],vMea,vRef);
    qEstE = KFEul(gyro,sf,q0,[],vMea,vRef);
    qEstQ = KFQua(gyro,sf,q0,[],vMea,vRef);
    qEstC = Comp(gyro,sf,q0,[],vMea,vRef);
    
    parsave(n,qTrue,gyro,vMea,vRef,qEstM,qEstE,qEstQ,qEstC,true);
end

parfor n = 1:N
    % true trajectory and measurement
    [qTrue,gyro] = genRotAxis(tTot,sf);
    vMea1 = genMea(qTrue,'v',vRef(:,1),'vstd',0.1);
    vMea2 = genMea(qTrue,'v',vRef(:,2),'vstd',0.1);
    vMea = cat(2,reshape(vMea1',3,1,[]),reshape(vMea2',3,1,[]));
    
    % filters
    q0 = mulQua(qTrue(1,:),expQua([pi,0,0]));
    qEstM = MEKF(gyro,sf,q0,[],vMea,vRef);
    qEstE = KFEul(gyro,sf,q0,[],vMea,vRef);
    qEstQ = KFQua(gyro,sf,q0,[],vMea,vRef);
    qEstC = Comp(gyro,sf,q0,[],vMea,vRef);
    
    parsave(n,qTrue,gyro,vMea,vRef,qEstM,qEstE,qEstQ,qEstC,false);
end

rmpath('..\rotation3d','..\Filters Without Bias',...
    '..\Generate Motions','..\Vector to Rotation');

end


function parsave(n,qTrue,gyro,vMea,vRef,qEstM,qEstE,qEstQ,qEstC,t)

if t
    save(strcat('C:\result-filterComp\vecMea_noBias\trig\',num2str(n),'.mat'),...
        'qTrue','gyro','vMea','vRef','qEstM','qEstE','qEstQ','qEstC');
else
    save(strcat('C:\result-filterComp\vecMea_noBias\rotAxis\',num2str(n),'.mat'),...
        'qTrue','gyro','vMea','vRef','qEstM','qEstE','qEstQ','qEstC');
end

end

