function [  ] = test_findKForComp(  )

addpath('..\rotation3d','..\Filters Without Bias','..\Generate Motions');

N = 60;
sf = 200;

% start from true attitude
Kt = [0.0005,0.001,0.0015,0.002,0.0025,0.003];

parfor n = 1:N
    [qTrue,gyro] = genTrig(60,sf);
    qMea = genMea(qTrue);
    q0 = qTrue(1,:);
    qEstM = MEKF(gyro,sf,q0,qMea,'P0',0);
    qEstC = zeros(length(gyro),4,length(Kt));
    for nk = 1:length(Kt)
        qEstC(:,:,nk) = Comp(gyro,sf,q0,qMea,'K',Kt(nk));
    end
    
    parsave(n,qTrue,gyro,qMea,qEstM,qEstC,true);
end

% start from false attitude
Kt = [0.01,0.02,0.03,0.04,0.05,0.06];

parfor n = 1:N
    [qTrue,gyro] = genTrig(60,sf);
    qMea = genMea(qTrue);
    q0 = mulQua(qTrue(1,:),expQua([pi,0,0]));
    qEstM = MEKF(gyro,sf,q0,qMea);
    qEstC = zeros(length(gyro),4,length(Kt));
    for nk = 1:length(Kt)
        qEstC(:,:,nk) = Comp(gyro,sf,q0,qMea,'K',Kt(nk));
    end
    
    parsave(n,qTrue,gyro,qMea,qEstM,qEstC,false);
end

rmpath('..\rotation3d','..\Filters Without Bias','..\Generate Motions');

end


function parsave(n,qTrue,gyro,qMea,qEstM,qEstC,t)

if t
    save(strcat('C:\result-filterComp\findKForComp\true\',num2str(n),'.mat'),...
            'qTrue','gyro','qMea','qEstM','qEstC');
else
    save(strcat('C:\result-filterComp\findKForComp\false\',num2str(n),'.mat'),...
            'qTrue','gyro','qMea','qEstM','qEstC');
end

end

