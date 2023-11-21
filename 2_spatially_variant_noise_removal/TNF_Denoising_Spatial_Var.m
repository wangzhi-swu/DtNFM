function [rI, Par]   =  TNF_Denoising_Spatial_Var( nI, I, Par )
rI  =  nI; 
[h, w, ch]  = size(rI);
Par.h=h;  Par.w=w;  Par.ch=ch;

Par = SearchNeighborIndex( Par );
NoiPat   =	Image2Patch( nI, Par ); 
Par.All =  size(NoiPat, 2); 
Sigma_arrCh = zeros(3*Par.ps^2, Par.All); 

for iter = 1 : Par.Iter
    rI_prev = rI;
    % iterative regularization 
    rI = rI + Par.delta * (nI - rI); 
    
    % image to patch
    CurPat = Image2Patch( rI, Par );

    % estimate local noise variance
    for c = 1:Par.ch
        TempSigma_arrCh = sqrt(max(0, repmat(Par.nSig(c)^2, 1, Par.All) - mean((NoiPat((c-1)*Par.ps2+1:c*Par.ps2, :) - CurPat((c-1)*Par.ps2+1:c*Par.ps2, :)).^2)));
        Sigma_arrCh((c-1)*Par.ps2+1:c*Par.ps2, :) = repmat(TempSigma_arrCh, [Par.ps2, 1]);
    end
    SigmaCol = sqrt(max(0, repmat(mean(Par.nSig.^2), 1, Par.All) - mean((NoiPat - CurPat).^2))); % do not know the peak
%     SigmaCol = sqrt(max(0, mean(Par.nSig0.^2) .* mean(Peak_Pat) - mean((NoiPat - CurPat).^2))); % know the peak
    
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
        if Par.image<10; imwrite(rI_prev./255, ['0',num2str(Par.image), '_', num2str(psnr_bst), '_c', num2str(Par.Constant), '_lmd', num2str(Par.lambda), '_rho', num2str(Par.rho), '_alp', num2str(Par.alpha), '_t', num2str(Par.t), '_it', num2str(iter), '_N', num2str(Par.nlsp), '.png']);
        else; imwrite(rI_prev./255, [num2str(Par.image), '_', num2str(psnr_bst), '_c', num2str(Par.Constant), '_lmd', num2str(Par.lambda), '_rho', num2str(Par.rho), '_alp', num2str(Par.alpha), '_t', num2str(Par.t), '_it', num2str(iter), '_N', num2str(Par.nlsp), '.png']);end
        break; 
    end
end
return;





