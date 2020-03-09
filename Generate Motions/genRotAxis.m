function [ q, gyro, bias ] = genRotAxis( t, sf, nobias, varargin )
% RA is rotation axis, RARA is another rotation axis that RA rotates about
% parameters: av: speed for the rotation about RA (rad/s)
%             avRA: speed for RA rotating about RARA (rad/s)
%             others: same as "getTrig"

time = (0:1/sf:t)';

% default and input parameters
E.initRA = [0;0;1];      % fixed
E.av = 6;
E.RARA = [0;1;0];        % fixed
E.avRA = 1;

U.GyroAngleRW = 0.1*pi/180;
U.GyroBiasRW = 0.01*pi/180;

if ~exist('nobias','var') || isempty(nobias)
    nobias = true;
end

[E,U] = parseVar(varargin,E,U);

% true state
R1 = @(t)[cos(t*E.avRA), 0, sin(t*E.avRA);
          0, 1, 0;
          -sin(t*E.avRA), 0, cos(t*E.avRA)];
R2 = @(t)[cos(t*E.av)/2-cos(2*t*E.avRA)/2+(cos(t*E.av)*cos(2*t*E.avRA))/2+1/2, -cos(t*E.avRA)*sin(t*E.av), -(sin(2*t*E.avRA)*(cos(t*E.av)-1))/2;
           cos(t*E.avRA)*sin(t*E.av), cos(t*E.av), -sin(t*E.avRA)*sin(t*E.av);
           -(sin(2*t*E.avRA)*(cos(t*E.av)-1))/2,  sin(t*E.avRA)*sin(t*E.av), cos(t*E.av)/2+cos(2*t*E.avRA)/2-(cos(t*E.av)*cos(2*t*E.avRA))/2+1/2];
R = @(t)R1(t)*R2(t);

dR11 = @(t)E.avRA*cos(t*E.av)*sin(t*E.avRA)-2*E.avRA*sin(t*E.avRA)+E.av*cos(t*E.avRA)*sin(t*E.av)+6*E.avRA*cos(t*E.avRA)^2*sin(t*E.avRA)-2*E.av*cos(t*E.avRA)^3*sin(t*E.av)-6*E.avRA*cos(t*E.avRA)^2*cos(t*E.av)*sin(t*E.avRA);
dR12 = @(t)2*E.avRA*sin(t*E.av)*sin(2*t*E.avRA)-E.av*cos(t*E.av)*cos(2*t*E.avRA);
dR13 = @(t)E.avRA*cos(t*E.avRA)*(cos(t*E.av)+2*cos(t*E.avRA)^2-2*cos(t*E.avRA)^2*cos(t*E.av))-sin(t*E.avRA)*(E.av*sin(t*E.av)+4*E.avRA*cos(t*E.avRA)*sin(t*E.avRA)-2*E.av*cos(t*E.avRA)^2*sin(t*E.av)-4*E.avRA*cos(t*E.avRA)*cos(t*E.av)*sin(t*E.avRA));
dR21 = @(t)E.av*cos(t*E.avRA)*cos(t*E.av)-E.avRA*sin(t*E.avRA)*sin(t*E.av);
dR22 = @(t)-E.av*sin(t*E.av);
dR23 = @(t)-E.avRA*cos(t*E.avRA)*sin(t*E.av)-E.av*cos(t*E.av)*sin(t*E.avRA);
dR31 = @(t)cos(t*E.avRA)*(4*E.avRA*cos(t*E.av)-5*E.avRA+6*E.avRA*cos(t*E.avRA)^2-6*E.avRA*cos(t*E.avRA)^2*cos(t*E.av)+2*E.av*cos(t*E.avRA)*sin(t*E.avRA)*sin(t*E.av));
dR32 = @(t)4*E.avRA*cos(t*E.avRA)^2*sin(t*E.av)-2*E.avRA*sin(t*E.av)+2*E.av*cos(t*E.avRA)*cos(t*E.av)*sin(t*E.avRA);
dR33 = @(t)E.avRA*sin(t*E.avRA)-2*E.avRA*cos(t*E.av)*sin(t*E.avRA)-2*E.av*cos(t*E.avRA)*sin(t*E.av)-6*E.avRA*cos(t*E.avRA)^2*sin(t*E.avRA)+2*E.av*cos(t*E.avRA)^3*sin(t*E.av)+6*E.avRA*cos(t*E.avRA)^2*cos(t*E.av)*sin(t*E.avRA);
dR = @(t)[dR11(t),dR12(t),dR13(t);dR21(t),dR22(t),dR23(t);dR31(t),dR32(t),dR33(t)];

q = zeros(length(time),4);
gyro = zeros(length(time),3);
for nt = 1:length(time)
    q(nt,:) = rot2qua(R(time(nt)))';
    gyro(nt,:) = vee(R(time(nt))'*dR(time(nt)));
end

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
    
    if strcmp(inputs{i},'av')
        E.av = inputs{i+1};
    elseif strcmp(inputs{i},'avRA')
        E.avRA = inputs{i+1};
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

