function [E_Img, Par]   =  TNF_Denoising_real( N_Img, O_Img, Par )
E_Img = N_Img;  % Estimated Image
[h, w, ch]  = size(E_Img);
Par.h = h;
Par.w = w;
Par.ch = ch;
Par = SearchNeighborIndex( Par );

NoiPat   =	Image2Patch( N_Img, Par );
Par.All =  size(NoiPat, 2);
Par.Sigma = sqrt(mean(Par.nSig0.^2));
Sigma = ones(Par.ch, length(Par.SelfIndex));
Innerloop = 3; bi_clock = true;
for iter = 1 : Par.Iter
    Par.iter = iter;
    rI_prev = E_Img;
    % iterative regularization
    E_Img =	E_Img + Par.delta * (N_Img - E_Img);
    CurPat = Image2Patch( E_Img, Par );

    if (iter==1)
        for c=1:3
            TempSigma_arrCh = Par.lambda * Par.nSig0(c) * Sigma(c, :);%
            Sigma_arrCh((c-1)*Par.ps2+1:c*Par.ps2, :) = repmat(TempSigma_arrCh, [Par.ps2, 1]);
        end
        SigmaCol = Par.lambda * sqrt(abs( repmat(Par.Sigma^2, 1, Par.All) - mean((NoiPat - CurPat).^2)));
    else
        for c=1:3
            TempSigma_arrCh = Par.lambda * Sigma(c, :);
            Sigma_arrCh((c-1)*Par.ps2+1:c*Par.ps2, :) = repmat(TempSigma_arrCh, [Par.ps2, 1]);
        end
        SigmaCol = Par.lambda *  SigmaPats;
    end
    Par.nlsp = Par.nlsp - 10*bi_clock;
    bi_clock = ~bi_clock;
    if mod(iter-1, Innerloop)==0
        NL_mat  =  Block_Matching_Real(CurPat, Par);
    end
    [Y_hat, W_hat, Sigma, SigmaPats]  =  TNF_Estimation_real( NL_mat, Sigma_arrCh, SigmaCol, CurPat, Par );
    SigmaPats = SigmaPats./(W_hat(1,:));
    SigmaPats(isnan(SigmaPats)) = 0;
    
    E_Img = PGs2Image(Y_hat, W_hat, Par);
    E_Img(E_Img>255) = 255;
    E_Img(E_Img<0.0) = 0.0;
    PSNR  =  csnr( O_Img, E_Img, 0, 0 );
    SSIM  =  cal_ssim( O_Img, E_Img, 0, 0 );
    fprintf( 'it %d, PSNR = %f  SSIM = %f\n', iter, PSNR, SSIM );
    Par.PSNR(iter, Par.image)  =  PSNR;
    
    if (iter>1 && Par.PSNR(iter-1, Par.image)>PSNR)
        psnr_bst = Par.PSNR(iter-1, Par.image);
        if Par.image<10; imwrite(rI_prev./255, ['0',num2str(Par.image), '_', num2str(psnr_bst), '_c', num2str(Par.Constant), '_lmd', num2str(Par.lambda), '_rho', num2str(Par.rho), '_alp', num2str(Par.alpha), '_t', num2str(Par.t), '_it', num2str(iter-1), '_N', num2str(Par.nlsp), '.png']);
        else; imwrite(rI_prev./255, [num2str(Par.image), '_', num2str(psnr_bst), '_c', num2str(Par.Constant), '_lmd', num2str(Par.lambda), '_rho', num2str(Par.rho), '_alp', num2str(Par.alpha), '_t', num2str(Par.t), '_it', num2str(iter-1), '_N', num2str(Par.nlsp), '.png']);end
        break; 
    end
end
return;