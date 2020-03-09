function [ qEst, P ] = MEKF( gyro, sf, q0, qMea, vMea, vRef, varargin )
% Attitude estimation using multiplicative extanded Kalman filter, see
% "Kalman Filtering for Spacecraft Attitude Estimation", Lefferts, 1982.
% Note: the quaternion convention is different from the above reference!
% input parameters: gyro: measured angular velocity
%                   sf: sampling frequency
%                   q0: initial quaternion
%                   qMea: attitude measurements
%                   vMea: vector measurements
%                   vRef: reference vector in global frame
% output parameters: P: covariance matrix for MEKF

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
    R = eye(3)*U.rvstd^2;
else
    [~,~,rvvar] = vMea2R(vMea(:,:,1),vRef,U.vMeaStd,'q');
    R = eye(3)*rvvar;
end

% pre-allocate memory
qEst = zeros(Nt,4);
P = zeros(3,3,Nt);
if strcmp(meaType,'V')
    qMea = zeros(Nt,4);
end

% initialize
qEst(1,:) = q0;
P(:,:,1) = eye(3)*U.P0^2;

% filter iteration
for nt = 2:Nt
    % integration
    av = 0.5*(gyro(nt-1,:)+gyro(nt,:));
    qEst(nt,:) = mulQua(qEst(nt-1,:),expQua(av*dt));
    
    % uncertainty propagation
    F = expRot(av*dt)';
    P(:,:,nt) = F*P(:,:,nt-1)*F'+Q;
    
    % if measurement comes in vector form, convert it first
    if strcmp(meaType,'V')
        [~,qMea(nt,:)] = vMea2R(vMea(:,:,nt),vRef,U.vMeaStd);
    end
    
    % update
    K = P(:,:,nt)*(P(:,:,nt)+R)^-1;
    dv = K*logQua(mulQua(invQua(qEst(nt,:)),qMea(nt,:)),'v')';
    qEst(nt,:) = mulQua(qEst(nt,:),expQua(dv)');
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

