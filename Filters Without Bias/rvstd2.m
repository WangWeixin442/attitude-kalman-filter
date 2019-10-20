function [ qstd, estd ] = rvstd2( rvstd )

Ns = 1000000;

rv = randn(Ns,3)*rvstd;

% fit a covariance for quaternion
q = expQua(rv);
qstd = cov(q);

% fit a covariance for Euler angles
e = qua2eul(q);
estd = cov(e);

end

