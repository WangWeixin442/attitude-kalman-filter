function [ ] = test_singleVec_noBias( )

addpath('..\rotation3d','..\Filters Without Bias','..\Generate Motions');

% global settings
N = 60;
sf = 200;
tTot = 300;
vRef = [0;0;1];

for n = 1:N
    % true trajectory and measurement
    [qTrue,gyro] = genTrig(tTot,sf,[],'magp',pi/6);
    vMea = genMea(qTrue,'v',vRef,'vstd',0.1);
    vMea = permute(vMea,[2,3,1]);
    
    % filters
    q0 = qTrue(1,:);
    qEst1 = MEKF(gyro,sf,q0,[],vMea,vRef,'errLocal',true);
    qEst2 = MEKF(gyro,sf,q0,[],vMea,vRef,'errLocal',false);
end

addpath('..\rotation3d','..\Filters Without Bias','..\Generate Motions');

end

