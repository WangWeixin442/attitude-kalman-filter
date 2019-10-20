function [  ] = test_MEKFvsEul_noBias(  )

addpath('..\rotation3d','..\Filters Without Bias','..\Generate Motions');

N = 60;
sf = 200;

eul = [0,pi/6,pi/3,pi/2];

for ne = 1:length(eul)
    parfor n = 1:N
        [qTrue,gyro] = genTrig(60,sf,[],'magp',eul(ne));
        qMea = genMea(qTrue);
        q0 = mulQua(qTrue(1,:),expQua([pi,0,0]));
        qEstM = MEKF(gyro,sf,q0,qMea);
        qEstE = KFEul(gyro,sf,q0,qMea);

        parsave(n,ne,qTrue,gyro,qMea,qEstM,qEstE);
    end
end

rmpath('..\Filters Without Bias','..\Generate Motions');

end


function parsave(n,ne,qTrue,gyro,qMea,qEstM,qEstE)

save(strcat('C:\result-filterComp\MEKFvsEul_noBias\',num2str(ne),'\',num2str(n),'.mat'),...
    'qTrue','gyro','qMea','qEstM','qEstE');

end

