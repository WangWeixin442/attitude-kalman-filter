function [ qMea ] = genMea( q, type, vref, varargin )

Nt = length(q);

% default paramters
U.rvstd = 0.1;

if ~exist('type','var') || isempty(type)
    type = 'A';
end
if strcmpi(type,'V') && (~exist('vref','var') || isempty(vref))
    
end

U = parseVar(varargin,U);

% generate random maesurement
if strcmpi(type,'A')
    dv = randn(Nt,3)*U.rvstd;
    qMea = mulQua(q,expQua(dv));
else
    error('vector measurement is under development');
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
    else
        error(strcat('No parameter specified by ',inputs{i}));
    end
    
    i = i+2;
end

end

