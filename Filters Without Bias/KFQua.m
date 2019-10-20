function [ qEst, P ] = KFQua( gyro, sf, q0, qMea, varargin )

dt = 1/sf;
Nt = length(gyro);

% default parameters
U.P0 = 100;
U.GyroAngleRW = 0.1*pi/180;
U.rvstd = 0.1;

U = parseVar(varargin,U);

% Q and R
Q = rvstd2(U.GyroAngleRW*sqrt(dt));
R = rvstd2(U.rvstd);

% pre-allocate memory
qEst = zeros(Nt,4);
P = zeros(4,4,Nt);

% initialize
qEst(1,:) = q0;
P(:,:,1) = rvstd2(U.P0);

% filter iteration
for nt = 2:Nt
    % integration
    av = 0.5*(gyro(nt-1,:)+gyro(nt,:));
    dq = expQua(av*dt);
    qEst(nt,:) = mulQua(qEst(nt-1,:),dq);
    
    % uncertainty propagation
    F = [dq(1),-dq(2),-dq(3),-dq(4);
         dq(2),dq(1),dq(4),-dq(3);
         dq(3),-dq(4),dq(1),dq(2);
         dq(4),dq(3),-dq(2),dq(1)];
    G = [qEst(nt,1),-qEst(nt,2),-qEst(nt,3),-qEst(nt,4)
         qEst(nt,2),qEst(nt,1),-qEst(nt,4),qEst(nt,3)
         qEst(nt,3),qEst(nt,4),qEst(nt,1),-qEst(nt,2)
         qEst(nt,4),-qEst(nt,3),qEst(nt,2),qEst(nt,1)];
    P(:,:,nt) = F*P(:,:,nt-1)*F'+G*Q*G';
    
    % update
    K = P(:,:,nt)*(P(:,:,nt)+R)^-1;
    if sqrt(sum((qEst(nt,:)-qMea(nt,:)).^2))>sqrt(2)
        qMea(nt,:) = -qMea(nt,:);
    end
    dq = qMea(nt,:)-qEst(nt,:);
    qEst(nt,:) = qEst(nt,:)+(K*dq')';
    qEst(nt,:) = qEst(nt,:)/sqrt(sum(qEst(nt,:).^2));
    P(:,:,nt) = (eye(4)-K)*P(:,:,nt);
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
    else
        error(strcat('No parameter specified by ',inputs{i}));
    end
    
    i = i+2;
end

end

