function [ mea ] = genMea( q, type, vref, varargin )
% For a given rotation motion, generate attitude or vector measurements
% input parameters: q: true quaternions
%                   type: 'A' for attitude; 'V' for vector
%                   vref: reference vector in global frame
% other parameters: rvstd: std for rotation vector (used when type='A')
%                   vstd: std for vector (used when type='V')

Nt = length(q);

% default paramters
U.rvstd = 0.1;
U.vstd = 0.5;

if ~exist('type','var') || isempty(type)
    type = 'A';
end
if strcmpi(type,'V') && (~exist('vref','var') || isempty(vref))
    vref = [0;0;9.8];
end

U = parseVar(varargin,U);

% generate random maesurement
if strcmpi(type,'A')
    dv = randn(Nt,3)*U.rvstd;
    mea = mulQua(q,expQua(dv));
elseif strcmpi(type,'V')
    RInv = invRot(qua2rot(q));
    noise = randn(Nt,3)*U.vstd;
    mea = zeros(Nt,3);
    for n = 1:Nt
        mea(n,:) = (RInv(:,:,n)*vref)'+noise(n,:);
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
    
    if strcmp(inputs{i},'rvstd')
        U.rvstd = inputs{i+1};
    elseif strcmp(inputs{i},'vstd')
        U.vstd = inputs{i+1};
    else
        error(strcat('No parameter specified by ',inputs{i}));
    end
    
    i = i+2;
end

end

