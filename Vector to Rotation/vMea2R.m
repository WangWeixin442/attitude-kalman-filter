function [ R, q, rvvar ] = vMea2R( vMea, vRef, vstd, method )
% convert vector measurements into attitude measurements
% see "Three-Axis Attitude Determination from Vector Observations", Shuster, 1981
% "Attitude Determination Using Vector Observations and the Singular Value Decomposition", Markley, 1988
% input parameters: vMea: vector measurements in IMU frame
%                   vRef: vector reference in global frame
%                   vstd: std of vector measurement
%                   method: 'q' or 'R'
% output parameters: rvvar: std for the error rotation vector
% vMea and vRef should be given as a 3-by-n matrix
% vstd shoulb be given as a 1-by-n matrix

% default method is rotation matrix
if ~exist('method','var') || isempty(method)
    method = 'R';
end

% F
sigmaTot = sum(vstd.^-2)^(-1/2);
w = sigmaTot^2./vstd.^2;
F = zeros(3);
for n = 1:size(vMea,2)
    F = F+w(n)*vRef(:,n)*vMea(:,n)';
end

if strcmpi(method,'q')
    % Quaternion method
    % K
    K = [trace(F),sum(w.*cross(vMea,vRef),2)';
        sum(w.*cross(vMea,vRef),2),F+F'-trace(F)*eye(3)];
    
    % qOpt
    [V,D] = eig(K);
    [~,ind] = max(diag(D));
    q = V(:,ind);
    
    % rvstd
    rvvar = zeros(3);
    for n = 1:size(vMea,2)
        rvvar = rvvar+w(n)*vRef(:,n)*vRef(:,n)';
    end
    rvvar = sigmaTot^2*(eye(3)-rvvar)^-1;
    
    % convert to rotation matrix
    R = qua2rot(q);
    
elseif strcmpi(method,'R')
    % rotation matrix method
    % proper singular value decomposition
    [U,S,V] = svd(F);
    S = S*diag([1,1,det(U*V)]);
    U = U*diag([1,1,det(U)]);
    V = V*diag([1,1,det(V)]);
    
    % Ropt
    R = U*V';
    
    % rvstd
    D = (trace(S)*eye(3)-S);
    rvvar = sigmaTot^2*U*(eye(3)-S)*D^-2*U';
    
    % convert to quaternion
    q = rot2qua(R);
end

end

