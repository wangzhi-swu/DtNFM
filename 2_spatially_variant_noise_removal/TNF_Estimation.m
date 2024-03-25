function [ Y_hat, W_hat] = TNF_Estimation( NL_mat, Sigma_arr, CurPat, Par )
    p23=Par.ps2ch;  All=Par.All;  key=size(NL_mat, 2); % 
    N = Par.nlsp;
    Y_hat = zeros([p23, All]);
    W_hat = zeros([p23, All]);
    for i   =   1 : key % For each keypatch group
        Y   =   CurPat(:, NL_mat(1:N, i)); 
        mY = repmat(mean( Y, 2 ), 1, N); % 
        Y   =   Y-mY;
        Sig_cs = Sigma_arr(:, NL_mat(1:N,i));
        
        [X] 	=   TNF_ADMM( Y, Sig_cs, Par); 
        
        Y_hat(:,NL_mat(1:N,i))  =  Y_hat(:,NL_mat(1:N,i)) + X + mY;
        W_hat(:,NL_mat(1:N,i))  =  W_hat(:,NL_mat(1:N,i)) + ones(p23, N);
    end
end