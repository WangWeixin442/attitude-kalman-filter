function [ qEst ] = Comp( gyro, sf, q0, qMea, varargin )

dt = 1/sf;
Nt = length(gyro);

% default
U.K = 0.03;
U = parseVar(varargin,U);

% proportional gain
K = U.K;

% pre-allocate memory
qEst = zeros(Nt,4);

% initialize
qEst(1,:) = q0;
dv = [0;0;0];

% filter iteration
for nt = 2:Nt
    % integration
    av = 0.5*(gyro(nt-1,:)+gyro(nt,:));
    qEst(nt,:) = mulQua(qEst(nt-1,:),expQua(av*dt+dv'));
    
    % update
    dv = K*logQua(mulQua(invQua(qEst(nt,:)),qMea(nt,:)),'v')';
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
    
    if strcmp(inputs{i},'K')
        U.K = inputs{i+1};
    else
        error(strcat('No parameter specified by ',inputs{i}));
    end
    
    i = i+2;
end

end

