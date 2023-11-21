function [rI_prev, Par, psnr_bst, ssim_bst]   =  TNF_Denoising( nI, I, Par )
rI  =  nI;   % Estimated Image
[h, w, ch]  = size(rI);
Par.h=h;  Par.w=w;  Par.ch=ch;

Par = SearchNeighborIndex( Par ); 

% noisy image to patch
NoiPat   =	Image2Patch( nI, Par ); 
Par.All =  size(NoiPat, 2); 
Sigma_arrCh = zeros(3*Par.ps^2, Par.All);
Par.Sigma = sqrt(mean(Par.nSig0.^2));
bi_clock = true;
for iter = 1 : Par.Iter
    Par.iter = iter;
    % iterative regularization 
    rI_prev = rI;
    rI = rI + Par.delta * (nI - rI);
    % image to patch
    CurPat = Image2Patch( rI, Par );
    % estimate local noise variance
    for c = 1:Par.ch
        TempSigma_arrCh = sqrt(max(0, repmat(Par.nSig0(c)^2, 1, Par.All) - mean((NoiPat((c-1)*Par.ps2+1:c*Par.ps2, :) - CurPat((c-1)*Par.ps2+1:c*Par.ps2, :)).^2)));
        Sigma_arrCh((c-1)*Par.ps2+1:c*Par.ps2, :) = repmat(TempSigma_arrCh, [Par.ps2, 1]);
    end
    SigmaCol = sqrt(max(0, repmat(Par.Sigma^2, 1, Par.All) - mean((NoiPat - CurPat).^2))); % Estimated Local Noise Level
    
    Par.nlsp = Par.nlsp - 10*bi_clock; 
    bi_clock = ~bi_clock;
    NL_mat  =  Block_Matching(CurPat, Par); 
    
    [Y_hat, W_hat]  =  TNF_Estimation( NL_mat, Sigma_arrCh, SigmaCol, CurPat, Par );
    rI = PGs2Image(Y_hat, W_hat, Par);
    rI(rI>255) = 255;
    rI(rI<0.0) = 0.0;
    
    PSNR   =  csnr( I, rI, 0, 0 );
    SSIM   =  cal_ssim( I, rI, 0, 0 );
    fprintf( 'Iter=%d, PSNR = %f  SSIM = %f\n', iter, PSNR, SSIM );
    Par.PSNR(iter, Par.image)  =  PSNR;
    Par.SSIM(iter, Par.image)  =  SSIM;

    if (iter>1 && Par.PSNR(iter-1, Par.image)>PSNR)
        psnr_bst = Par.PSNR(iter-1, Par.image);
        ssim_bst = Par.SSIM(iter-1, Par.image);
        if Par.image<10; imwrite(rI_prev./255, ['0',num2str(Par.image), '_', num2str(psnr_bst), '_lmd', num2str(Par.lambda), '_rho', num2str(Par.rho), '_alp', num2str(Par.alpha), '_t', num2str(Par.t), '_it', num2str(iter), '_N', num2str(Par.nlsp), '.png']);
        else; imwrite(rI_prev./255, [num2str(Par.image), '_', num2str(psnr_bst), '_lmd', num2str(Par.lambda), '_rho', num2str(Par.rho), '_alp', num2str(Par.alpha), '_t', num2str(Par.t), '_it', num2str(iter), '_N', num2str(Par.nlsp), '.png']);end
        break; 
    end
end
return;





