function [rI, psnr_f, ssim_f]   =  TNF_Denoising( nI, I, Par )
rI  =  nI;   % Estimated Image
[h, w, ch]  = size(rI);
Par.h=h;  Par.w=w;  Par.ch=ch;

Par = SearchNeighborIndex( Par ); 

% noisy image to patch
NoiPat   =	Image2Patch( nI, Par ); 
Par.All =  size(NoiPat, 2); 
Sigma_arrCh = zeros(3*Par.ps^2, Par.All); % 3p^2 * All
Par.Sigma = sqrt(mean(Par.nSig0.^2));
bi_clock = true;
for iter = 1 : Par.Iter
    Par.iter = iter;
    
    rI_prev = rI;
    rI = rI + Par.delta * (nI - rI); % iterative regularization 
    % image to patch
    CurPat = Image2Patch( rI, Par );
    % estimate local noise variance
    for c = 1:Par.ch
        TempSigma_arrCh = sqrt(abs( repmat(Par.nSig0(c)^2, 1, Par.All) - mean((NoiPat((c-1)*Par.ps2+1:c*Par.ps2, :) - CurPat((c-1)*Par.ps2+1:c*Par.ps2, :)).^2)));
        Sigma_arrCh((c-1)*Par.ps2+1:c*Par.ps2, :) = Par.lambda .* repmat(TempSigma_arrCh, [Par.ps2, 1]);
    end
    if iter==1
        Sigma_arrCh = Sigma_arrCh ./ Par.lambda;
    end
    
    Par.nlsp = Par.nlsp - 10*bi_clock; % 60,60, 50,50, 40, ...
    bi_clock = ~bi_clock;
    NL_mat  =  Block_Matching(CurPat, Par); % Caculate Non-local similar patches for each
    
    [Y_hat, W_hat]  =  TNF_Estimation( NL_mat, Sigma_arrCh, CurPat, Par );
    rI = PGs2Image(Y_hat, W_hat, Par);
    rI(rI>255) = 255;
    rI(rI<0.0) = 0.0;
    
    PSNR   =  csnr( I, rI, 0, 0 );
    SSIM   =  cal_ssim( I, rI, 0, 0 );
    fprintf( 'Iter=%d, PSNR = %f  SSIM = %f\n', iter, PSNR, SSIM );
    Par.PSNR(iter, Par.image)  =  PSNR;
    Par.SSIM(iter, Par.image)  =  SSIM;

    if (iter>1 && Par.PSNR(iter-1, Par.image)>PSNR)
        rI = rI_prev;
        psnr_f = Par.PSNR(iter-1, Par.image);
        ssim_f = Par.SSIM(iter-1, Par.image);
        if Par.image<10; imwrite(rI./255, ['0',num2str(Par.image), '_', num2str(psnr_f), '_lmd', num2str(Par.lambda), '_rho', num2str(Par.rho), '_alp', num2str(Par.alpha), '_t', num2str(Par.t), '_it', num2str(iter), '_N', num2str(Par.nlsp), '.png']);
        else; imwrite(rI./255, [num2str(Par.image), '_', num2str(psnr_f), '_lmd', num2str(Par.lambda), '_rho', num2str(Par.rho), '_alp', num2str(Par.alpha), '_t', num2str(Par.t), '_it', num2str(iter), '_N', num2str(Par.nlsp), '.png']);end
        break; 
    end
    if iter==Par.Iter
        psnr_f = PSNR; ssim_f = SSIM;
        if Par.image<10; imwrite(rI./255, ['0',num2str(Par.image), '_', num2str(PSNR), '_lmd', num2str(Par.lambda), '_rho', num2str(Par.rho), '_alp', num2str(Par.alpha), '_t', num2str(Par.t), '_it', num2str(iter), '_N', num2str(Par.nlsp), '.png']);
        else; imwrite(rI./255, [num2str(Par.image), '_', num2str(PSNR), '_lmd', num2str(Par.lambda), '_rho', num2str(Par.rho), '_alp', num2str(Par.alpha), '_t', num2str(Par.t), '_it', num2str(iter), '_N', num2str(Par.nlsp), '.png']);end
        break; 
    end
end
return;





