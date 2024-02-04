function [SigmaZ, svp] = Closed_TL12 (SigmaY, lambda, alpha, C, t)
    Lambda = lambda*sqrt(C); % Lambda ∝ lambda/rho
    n = length(SigmaY);
    SigmaZ = zeros(n, 1);

    z = max(SigmaY(t+1:n) - Lambda, 0); % z \in R^{n-t}
    norm_2 = sqrt( sum(z.^2) )+eps;
    K = 1 + alpha*Lambda/norm_2;

    SigmaZ(1:t) = SigmaY(1:t);
    v = 0; % this makes 1) the threshold Lambda -> Lambda - v/K; 
    % and 2) the output singular values (shrunk by the new Lambda) where 0 <= SigmaZ <= v be imposed to 0. 
    % This change will not break the convergence of the algorithm. 
    SigmaZ(t+1:n) = K*z + v;

    over = find( SigmaY(t+1:n) > (K*Lambda - v )/(K-1) ) + t;
    SigmaZ(over) = SigmaY(over);

    svp = length(find(z > 0)) + t;
    SigmaZ = SigmaZ(1:svp);
end
% plot(SigmaY,SigmaY, 'k*-');hold on; plot(SigmaY(1:svp),SigmaZ(1:svp),'bs--');hold off;
