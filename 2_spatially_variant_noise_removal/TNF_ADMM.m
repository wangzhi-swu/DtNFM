function  [X] =  TNF_ADMM( Y, NSig, SigCol, Par )
p23 = Par.ps2ch; N = Par.nlsp;
X = zeros([p23, N]);
Z = zeros([p23, N]);
A = zeros([p23, N]);

%% 双加权
var1=(std(NSig))/(mean(NSig) + eps);
var2=(std(SigCol))/(mean(SigCol) + eps);
p = (var1+eps)/(var1+var2 +2*eps);
%% 
mNSig = min( NSig([1,floor(end/2),end]) ) +eps;
mScol = min( SigCol ) +eps;
W = mNSig^(p) ./ (NSig.^(p) + eps); % 
v = mScol^(1-p) ./ (SigCol.^(1-p) + eps); %

rho = Par.rho;
mu = Par.mu;
lambda = Par.lambda;
alpha = Par.alpha;
c = Par.Constant;

for iter = 1 : Par.maxIter
    % update X
    X = Y + (rho.*(Z-Y) -A) ./ (2.*W.^2 * v.^2 + rho); % equivalent to X = (Y.*(W.^2*v.^2) + 0.5*rho*Z - 0.5*A) ./ (W.^2*v.^2 + 0.5*rho);

    % update Z
    [U, SigmaY, V] = svd(full(X + A/rho), 'econ');
    C = c*sqrt(N)*(mScol^4 + mNSig^4)/(rho^2); % spatially variant noise
    [SigmaZ , svp] = Closed_TL12(diag(SigmaY), lambda, alpha, C, Par.t);
    Z = U(:, 1:svp) * diag(SigmaZ) * V(:, 1:svp)';

    % update A and rho
    A = A + rho*(X - Z);
    rho = min(1e4, mu*rho);
end
return;