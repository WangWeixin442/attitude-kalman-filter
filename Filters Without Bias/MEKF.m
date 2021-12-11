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
    num_vec = size(vMea,2);
else
    meaType = 'A';
end

% default parameters
U.P0 = 100;
U.GyroAngleRW = 0.1*pi/180;
U.rvstd = 0.1;
U.errLocal = true;

if strcmpi(meaType,'V')
    U.vMeaStd = [0.1,0.1];
    U.vec2R = true;
end

U = parseVar(varargin,U);

% Q and R
Q = eye(3)*U.GyroAngleRW^2*dt;
if strcmp(meaType,'A')
    R = eye(3)*U.rvstd^2;
elseif num_vec>1 && U.vec2R
    [~,~,rvvar] = vMea2R(vMea(:,:,1),vRef,U.vMeaStd,'q');
    R = eye(3)*rvvar;
else
    R = zeros(3*num_vec,3*num_vec);
    for i = 1:num_vec
        R(3*(i-1)+1:3*i,3*(i-1)+1:3*i) = U.vMeaStd(i)^2*eye(3);
    end
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
    if U.errLocal
        F = expRot(av*dt)';
    else
        F = eye(3);
    end
    P(:,:,nt) = F*P(:,:,nt-1)*F'+Q;
    
    % measurement update
    if strcmp(meaType,'V') && (num_vec==1 || ~vec2R)
        % if measurement comes in vector form, and if there is only one
        % vector measurement or when vec2R is false, we directly use
        % linearization of the measurement function to update
        H = zeros(3*num_vec,3);
        for i = 1:num_vec
            if U.errLocal
                H(3*(i-1)+1:3*i,:) = hat(qua2rot(qEst(nt,:))'*vRef(:,i));
            else
                H(3*(i-1)+1:3*i,:) = qua2rot(qEst(nt,:))'*hat(vRef(:,i));
            end
        end
        
        K = P(:,:,nt)*H'*(H*P(:,:,nt)*H'+R)^-1;
        dv = K*reshape(vMea(:,:,nt)-qua2rot(qEst(nt,:))'*vRef,3*num_vec,1);
        if U.errLocal
            qEst(nt,:) = mulQua(qEst(nt,:),expQua(dv)');
        else
            qEst(nt,:) = mulQua(expQua(dv)',qEst(nt,:));
        end
        P(:,:,nt) = (eye(3)-K*H)*P(:,:,nt);
    else
        if strcmp(meaType,'V')
            % if measurement comes in vector form, convert it first
            [~,qMea(nt,:)] = vMea2R(vMea(:,:,nt),vRef,U.vMeaStd);
        end
        
        K = P(:,:,nt)*(P(:,:,nt)+R)^-1;
        if U.errLocal
            dv = K*logQua(mulQua(invQua(qEst(nt,:)),qMea(nt,:)),'v')';
            qEst(nt,:) = mulQua(qEst(nt,:),expQua(dv)');
        else
            dv = K*logQua(mulQua(qMea(nt,:),invQua(qEst(nt,:))),'v')';
            qEst(nt,:) = mulQua(expQua(dv)',qEst(nt,:));
        end
        P(:,:,nt) = (eye(3)-K)*P(:,:,nt);
    end
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
    
    if strcmpi(inputs{i},'P0')
        U.P0 = inputs{i+1};
    elseif strcmpi(inputs{i},'GyroAngleRW')
        U.GyroAngleRW = inputs{i+1};
    elseif strcmpi(inputs{i},'rvstd')
        U.rvstd = inputs{i+1};
    elseif strcmpi(inputs{i},'vMeaStd')
        U.vMeaStd = inputs{i+1};
    elseif strcmpi(inputs{i},'errLocal')
        U.errLocal = inputs{i+1};
    elseif strcmpi(inputs{i},'vec2R')
        U.vec2R = inputs{i+1};
    else
        error(strcat('No parameter specified by ',inputs{i}));
    end
    
    i = i+2;
end

end

