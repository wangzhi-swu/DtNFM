function  [Z, sigma, sigPats] =  TNF_ADMM_real( Y, NSig, SigCol, Par )
mNSig = min(NSig([1, floor(end/2), end]));
mScol = min(SigCol);
mean_NSig = mean(NSig([1, floor(end/2), end]));
mean_SigCol = mean(SigCol);
var1 = std(NSig) / (mean_NSig + eps);
var2 = std(SigCol) / (mean_SigCol + eps);
p=(var1+eps)/(var1+var2 +2*eps);
p = max(p,0.95);

W = mNSig^p ./ (NSig.^p + eps);
v = mScol^(1-p) ./ (SigCol.^(1-p) + eps);
rho = Par.rho;
mu = Par.mu;
lambda = Par.lambda;
alpha = Par.alpha;
c = Par.Constant;
N = Par.nlsp;

Z = zeros(size(Y));
A = zeros(size(Y));
iter = 0;
while iter < Par.maxIter
    iter = iter + 1;
    % update X, fix Z and A
    X = Y + (rho.*(Z-Y) -A) ./ (2.*W.^2 * v.^2 + rho); % equivalent to  X = (Y.*(W.^2*v.^2) + 0.5*rho*Z - 0.5*A) ./ (W.^2*v.^2 + 0.5*rho);
    
    % update Z, fix X and A
    X_A = X + A/rho;
    [U, SigmaY, V] = svd(full(X_A), 'econ');
    C = c*sqrt(N)*( mNSig^2 ) / rho^2;%
    [SigmaZ , svp] = Closed_TL12(diag(SigmaY), lambda, alpha, C, Par.t);
    Z =  U(:, 1:svp) * diag(SigmaZ) * V(:, 1:svp)';
    % update the multiplier A, fix Z and X
    A = A + rho * (X - Z);
    rho = mu * rho;
end
sigma = sqrt( mean( reshape(mean((Y - Z).^2, 2), [Par.ps2 Par.ch])) )';
sigPats = sqrt( mean( (Y - Z).^2 ) );
return;
