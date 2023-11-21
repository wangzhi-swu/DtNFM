function [ Y_hat, W_hat, Sigma, SigmaPats] = TNF_Estimation_real( NL_mat, Sigma_arr, Sigma_Col, CurPat, Par )
Y_hat = zeros(size(CurPat));
W_hat = zeros(size(CurPat));
Sigma = zeros(Par.ch, length(Par.SelfIndex));
SigmaPats = zeros(1, Par.All); 
for i   =   1 : length(Par.SelfIndex) % For each keypatch group
    Y   =   CurPat(:, NL_mat(1:Par.nlsp,i)); % Non-local similar patches to the keypatch
    mY  =   repmat(mean( Y, 2 ),1,Par.nlsp);
    Y   =   Y-mY;
    
    [X, Sigma(:, i), sigPats] 	=   TNF_ADMM_real( Y, Sigma_arr(:, i), Sigma_Col(:, NL_mat(1:Par.nlsp,i)), Par); 
    
    Y_hat(:,NL_mat(1:Par.nlsp,i))  =  Y_hat(:,NL_mat(1:Par.nlsp,i))+X+mY;
    W_hat(:,NL_mat(1:Par.nlsp,i))  =  W_hat(:,NL_mat(1:Par.nlsp,i))+ones(Par.ps2ch, Par.nlsp);
    SigmaPats(NL_mat(1:Par.nlsp,i)) = SigmaPats(NL_mat(1:Par.nlsp,i)) + sigPats;
end
end