function [ qEst, P ] = KFEul( gyro, sf, q0, qMea, vMea, vRef, varargin )

dt = 1/sf;
Nt = length(gyro);

% rotation or vector measurement
if isempty(qMea)
    meaType = 'V';
else
    meaType = 'A';
end

% default parameters
U.P0 = 100;
U.GyroAngleRW = 0.1*pi/180;
U.rvstd = 0.1;
U.vMeaStd = [0.1,0.1];

U = parseVar(varargin,U);

% Q and R
Q = eye(3)*U.GyroAngleRW^2*dt;
if strcmp(meaType,'A')
    [~,R] = rvstd2(U.rvstd);
else
    [~,~,rvvar] = vMea2R(vMea(:,:,1),vRef,U.vMeaStd,'q');
    [~,R] = rvstd2(sqrt(rvvar));
end

% pre-allocate memory
qEst = zeros(Nt,4);
P = zeros(3,3,Nt);
if strcmp(meaType,'V')
    qMea = zeros(Nt,4);
end

% initialize
qEst(1,:) = q0;
[~,P(:,:,1)] = rvstd2(U.P0);

% filter iteration
for nt = 2:Nt
    % integration
    av = 0.5*(gyro(nt-1,:)+gyro(nt,:));
    dq = expQua(av*dt);
    qEst(nt,:) = mulQua(qEst(nt-1,:),dq);
    
    % uncertainty propagation
    e = qua2eul(qEst(nt,:));
    WB(1,1:3) = [1, sin(e(1))*tan(e(2)), cos(e(1))*tan(e(2))];
    WB(2,1:3) = [0, cos(e(1)), -sin(e(1))];
    WB(3,1:3) = [0, sin(e(1))/cos(e(2)), cos(e(1))/cos(e(2))];
    VB(1,1) = cos(e(1))*tan(e(2))*gyro(2) - sin(e(1))*tan(e(2))*gyro(3);
    VB(1,2) = sin(e(1))*cos(e(2))^-2*gyro(2) + cos(e(1))*cos(e(2))^-2*gyro(3);
    VB(1,3) = 0;
    VB(2,1) = -sin(e(1))*gyro(2)-cos(e(1))*gyro(3);
    VB(2,2) = 0;
    VB(2,3) = 0;
    VB(3,1) = cos(e(1))/cos(e(2))*gyro(2) - sin(e(1))/cos(e(2))*gyro(3);
    VB(3,2) = sin(e(1))*sin(e(2))*cos(e(2))^-2*gyro(2) + cos(e(1))*sin(e(2))*cos(e(2))^-2*gyro(3);
    VB(3,3) = 0;
    
    F = eye(3)+VB*dt;
    P(:,:,nt) = F*P(:,:,nt-1)*F'+WB*Q*WB';
    
    % if measurement comes in vector form, convert it first
    if strcmp(meaType,'V')
        [~,qMea(nt,:)] = vMea2R(vMea(:,:,nt),vRef,U.vMeaStd);
    end
    
    % update
    K = P(:,:,nt)*(P(:,:,nt)+R)^-1;
    de = K*wrapToPi(qua2eul(qMea(nt,:))-e)';
    qEst(nt,:) = eul2qua(e+de');
    P(:,:,nt) = (eye(3)-K)*P(:,:,nt);
end

end


function [ U ] = parseVar( inputs, U )

i = 1;
while i <= length(inputs)
    % not name-value pairs
    if i == length(inputs)
        error(strcat('No value assigned to ',inputs{i}));
    end
    
    % first argument is not string
    if ~ischar(inputs{i})
        error('Name should be a string');
    end
    
    if strcmp(inputs{i},'P0')
        U.P0 = inputs{i+1};
    elseif strcmp(inputs{i},'GyroAngleRW')
        U.GyroAngleRW = inputs{i+1};
    elseif strcmp(inputs{i},'rvstd')
        U.rvstd = inputs{i+1};
    elseif strcmp(inputs{i},'vMeaStd')
        U.vMeaStd = inputs{i+1};
    else
        error(strcat('No parameter specified by ',inputs{i}));
    end
    
    i = i+2;
end

end

