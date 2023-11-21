clear; 

GT_image_dir  =    'kodim_spatial_var/'; 
Noi_image_dir  =    'SV_Noi_35/'; 
GT_fpath = fullfile(GT_image_dir, '*.png');
Noi_fpath = fullfile(Noi_image_dir, '*.png');
GT_im_dir  = dir(GT_fpath);
Noi_im_dir = dir(Noi_fpath);
im_num = length(GT_im_dir);

Sig_Level = [30 35 40];

for Truncate = [8]  
for ALPHA = [1]
for LAMBDA = [1.5]
for RHO = [1.6]
% parameters for denoising
Par.nSig0     =   Sig_Level; % STD of the noise image
Par.win       =   20;        % Non-local patch searching window
Par.Constant  =   1;         % Constant num for the weight vector
Par.Innerloop =   2;         % InnerLoop Num of between re-blockmatching
Par.ps        =   5;         % Patch size, larger values would get better performance, but will be slower
Par.step      =   4;         % The step between neighbor patches, smaller values would get better performance, but will be slower
Par.Iter      =   12;        % iteration upperbound
Par.delta     =   0.1;       % step back length in iterative regularization
Par.maxIter   =   10;

Par.rho       =   RHO; 
Par.mu        =   1.002;

% this parameter is not finally determined yet
Par.lambda    =   LAMBDA;%0.85;          % for different noise levels, this parameter should be tuned to achieve better performance
Par.t = Truncate; 


% record all the results in each iteration
Par.PSNR = zeros(Par.Iter, im_num, 'single');
Par.SSIM = zeros(Par.Iter, im_num, 'single');
PSNRs = zeros(im_num,1); SSIMs = zeros(im_num,1);
for i = [10 15 20 24] 
    Par.alpha = ALPHA;
    Par.image   =   i;
    Par.nlsp    =   60; 
    Par.I = double( imread(fullfile(GT_image_dir, GT_im_dir(i).name)) );
    S = regexp(GT_im_dir(i).name, '\.', 'split');
    [h, w, ch] = size(Par.I);
    
    noise_peak = peaks(h);
    M = max(noise_peak(:));
    m = min(noise_peak(:));
    noise_peak = (noise_peak-m)./(M-m);
    noise_peak = flipud(noise_peak);
    Par.peak = noise_peak; 
    Par.nSig = Sig_Level .* mean( noise_peak(:) );
    Par.nim = double( imread(fullfile(Noi_image_dir, Noi_im_dir(i).name)) );
    fprintf('%s c=%.2f lambda=%.2f rho=%.2f alpha=%.2f t=%d N=%d\n', GT_im_dir(i).name(6:7), Par.Constant, Par.lambda, Par.rho, Par.alpha, Par.t, Par.nlsp);
    PSNR  =   csnr( Par.nim, Par.I, 0, 0 );
    SSIM  =   cal_ssim( Par.nim, Par.I, 0, 0 );
    fprintf('  PSNR_0 = %f  SSIM_0 = %f\n', PSNR, SSIM);
    PSNRs(i) = PSNR; SSIMs(i) = SSIM;
    time0 = clock;
    [im_out, Par] = TNF_Denoising_Spatial_Var( Par.nim, Par.I, Par );
    fprintf('Time used %f s\n',(etime(clock, time0)));
end
end
end
fprintf('\n\n');
end
end
