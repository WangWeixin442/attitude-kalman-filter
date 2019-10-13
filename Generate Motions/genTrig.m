function [ q, gyro, bias ] = genTrig( t, sf, nobias, varargin )

time = (0:1/sf:t)';

% default and input parameters
E.fr = 0.35; E.fp = 0.35; E.fh = 0.35;
E.phir = 0; E.phip = 0; E.phih = 0;
E.magr = pi; E.magp = pi/2; E.magh = pi;
E.medr = 0; E.medp = 0; E.medh = 0;

U.GyroAngleRW = 0.1*pi/180;
U.GyroBiasRW = 0.01*pi/180;

if ~exist('nobias','var') || isempty(nobias)
    nobias = true;
end

[E,U] = parseVar(varargin,E,U);

% true state
roll = @(t)E.magr*sin(E.fr*2*pi*t+E.phir)+E.medr;
pitch = @(t)E.magp*sin(E.fp*2*pi*t+E.phip)+E.medp;
yaw = @(t)E.magh*sin(E.fh*2*pi*t+E.phih)+E.medh;

sr2 = @(t)sin(roll(t)/2); sp2 = @(t)sin(pitch(t)/2); sy2 = @(t)sin(yaw(t)/2);
cr2 = @(t)cos(roll(t)/2); cp2 = @(t)cos(pitch(t)/2); cy2 = @(t)cos(yaw(t)/2);
q1 = @(t)cr2(t).*cp2(t).*cy2(t)+sr2(t).*sp2(t).*sy2(t);
q2 = @(t)sr2(t).*cp2(t).*cy2(t)-cr2(t).*sp2(t).*sy2(t);
q3 = @(t)cr2(t).*sp2(t).*cy2(t)+sr2(t).*cp2(t).*sy2(t);
q4 = @(t)cr2(t).*cp2(t).*sy2(t)-sr2(t).*sp2(t).*cy2(t);
q = [q1(time),q2(time),q3(time),q4(time)];

dr = @(t)E.magr*E.fr*2*pi*cos(E.fr*2*pi*t+E.phir);
dp = @(t)E.magp*E.fp*2*pi*cos(E.fp*2*pi*t+E.phip);
dy = @(t)E.magh*E.fh*2*pi*cos(E.fh*2*pi*t+E.phih);

dq1 = @(t)-sr2(t).*cp2(t).*cy2(t)*0.5.*dr(t)-cr2(t).*sp2(t).*cy2(t)*0.5.*dp(t)-cr2(t).*cp2(t).*sy2(t)*0.5.*dy(t)...
          +cr2(t).*sp2(t).*sy2(t)*0.5.*dr(t)+sr2(t).*cp2(t).*sy2(t)*0.5.*dp(t)+sr2(t).*sp2(t).*cy2(t)*0.5.*dy(t);
dq2 = @(t)+cr2(t).*cp2(t).*cy2(t)*0.5.*dr(t)-sr2(t).*sp2(t).*cy2(t)*0.5.*dp(t)-sr2(t).*cp2(t).*sy2(t)*0.5.*dy(t)...
          +sr2(t).*sp2(t).*sy2(t)*0.5.*dr(t)-cr2(t).*cp2(t).*sy2(t)*0.5.*dp(t)-cr2(t).*sp2(t).*cy2(t)*0.5.*dy(t);
dq3 = @(t)-sr2(t).*sp2(t).*cy2(t)*0.5.*dr(t)+cr2(t).*cp2(t).*cy2(t)*0.5.*dp(t)-cr2(t).*sp2(t).*sy2(t)*0.5.*dy(t)...
          +cr2(t).*cp2(t).*sy2(t)*0.5.*dr(t)-sr2(t).*sp2(t).*sy2(t)*0.5.*dp(t)+sr2(t).*cp2(t).*cy2(t)*0.5.*dy(t);
dq4 = @(t)-sr2(t).*cp2(t).*sy2(t)*0.5.*dr(t)-cr2(t).*sp2(t).*sy2(t)*0.5.*dp(t)+cr2(t).*cp2(t).*cy2(t)*0.5.*dy(t)...
          -cr2(t).*sp2(t).*cy2(t)*0.5.*dr(t)-sr2(t).*cp2(t).*cy2(t)*0.5.*dp(t)+sr2(t).*sp2(t).*sy2(t)*0.5.*dy(t);

w1 = @(t)(q1(t).*dq1(t)+q2(t).*dq2(t)+q3(t).*dq3(t)+q4(t).*dq4(t))*2;
w2 = @(t)(q1(t).*dq2(t)-q2(t).*dq1(t)-q3(t).*dq4(t)+q4(t).*dq3(t))*2;
w3 = @(t)(q1(t).*dq3(t)+q2(t).*dq4(t)-q3(t).*dq1(t)-q4(t).*dq2(t))*2;
w4 = @(t)(q1(t).*dq4(t)-q2(t).*dq3(t)+q3(t).*dq2(t)-q4(t).*dq1(t))*2;

gyro = [w2(time),w3(time),w4(time)];

% corruption of inertial errors
if ~nobias
    bias = [zeros(1,3);cumsum(randn(length(time)-1,3)*U.GyroBiasRW/sqrt(sf))];
else
    bias = zeros(length(time),3);
end

gyroRW = randn(length(time),3)*U.GyroAngleRW*sqrt(sf);
gyro = gyro+gyroRW+bias;

end

function [ E, U ] = parseVar( inputs, E, U )

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
    
    if strcmp(inputs{i},'fgyro')
        E.fr = inputs{i+1}; E.fp = inputs{i+1}; E.fh = inputs{i+1};
    elseif strcmp(inputs{i},'fr')
        E.fr = inputs{i+1};
    elseif strcmp(inputs{i},'fp')
        E.fp = inputs{i+1};
    elseif strcmp(inputs{i},'fh')
        E.fh = inputs{i+1};
    elseif strcmp(inputs{i},'phigyro')
        E.phir = inputs{i+1}; E.phip = inputs{i+1}; E.phih = inputs{i+1};
    elseif strcmp(inputs{i},'phir')
        E.phir = inputs{i+1};
    elseif strcmp(inputs{i},'phip')
        E.phip = inputs{i+1};
    elseif strcmp(inputs{i},'phih')
        E.phih = inputs{i+1};
    elseif strcmp(inputs{i},'magr')
        if abs(inputs{i+1})>pi
            error('magnitude of roll must be within [0,pi]')
        end
        E.magr = inputs{i+1};
    elseif strcmp(inputs{i},'magp')
        if abs(inputs{i+1})>pi
            error('magnitude of pitch must be within [0,pi/2]')
        end
        E.magp = inputs{i+1};
    elseif strcmp(inputs{i},'magh')
        if abs(inputs{i+1})>pi
            error('magnitude of yaw must be within [0,pi]')
        end
        E.magh = inputs{i+1};
    elseif strcmp(inputs{i},'medr')
        if inputs{i+1}+E.magr>pi || inputs{i+1}-E.magr<-pi
            error('roll angle must be within [-pi,pi]')
        end
        E.medr = inputs{i+1};
    elseif strcmp(inputs{i},'medp')
        if inputs{i+1}+E.magp>pi/2 || inputs{i+1}-E.magp<-pi/2
            error('pitch angle must be within [-pi/2,pi/2]')
        end
        E.medp = inputs{i+1};
    elseif strcmp(inputs{i},'medh')
        if inputs{i+1}+E.magy>pi || inputs{i+1}-E.magy<-pi
            error('yaw angle must be within [-pi,pi]')
        end
        E.medh = inputs{i+1};
    elseif strcmp(inputs{i},'GyroAngleRW')
        U.GyroAngleRW = inputs{i+1};
    elseif strcmp(inputs{i},'GyroBiasRW')
        U.GyroBiasRW = inputs{i+1};
    else
        error(strcat('No parameter specified by ',inputs{i}));
    end
    
    i = i+2;
end

end

