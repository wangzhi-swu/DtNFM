clear; 

GT_image_dir  =    'kodim_spatial_var/'; 
Noi_image_dir  =    'SV_Noi_35/'; 
GT_fpath = fullfile(GT_image_dir, '*.png');
Noi_fpath = fullfile(Noi_image_dir, '*.png');
GT_im_dir  = dir(GT_fpath);
Noi_im_dir = dir(Noi_fpath);
im_num = length(GT_im_dir);

Sig_Level = [30 35 40];
for Truncate = [2]  
for ALPHA = [1.5]
for LAMBDA = [ 0.8]
%     lb = -LAMBDA+3.1;
for RHO = [0.45]
% if RHO<=0.16; continue; end
% parameters for denoising
Par.nSig0     =   Sig_Level; % STD of the noise image
Par.win       =   20;        % Non-local patch searching window
Par.Constant  =   sqrt(2);   % Constant num for the weight vector
Par.Innerloop =   2;         % InnerLoop Num of between re-blockmatching
Par.ps        =   4;         % Patch size, larger values would get better performance, but will be slower
Par.step      =   3;         % The step between neighbor patches, smaller values would get better performance, but will be slower
Par.Iter      =   12;        % total iter numbers
Par.delta     =   0.1;
Par.maxIter   =   10;
Par.rho       =   RHO;       % In final version, this parameter will be changed
Par.mu        =   1.002;
Par.lambda    =   LAMBDA;    % for different noise levels, this parameter should be tuned to achieve better performance
Par.t         =   Truncate; 

% record all the results in each iteration
Par.PSNR = zeros(Par.Iter, im_num, 'single');
Par.SSIM = zeros(Par.Iter, im_num, 'single');
PSNRs = zeros(im_num,1); SSIMs = zeros(im_num,1); Times = zeros([im_num,1]);
for i = 1
    Par.alpha = ALPHA;
    Par.image   =   i;
    Par.nlsp    =   70;   % Initial Non-local Similar Patches number
    Par.I = double( imread(fullfile(GT_image_dir, GT_im_dir(i).name)) );
    S = regexp(GT_im_dir(i).name, '\.', 'split');
    [h, w, ch] = size(Par.I);
    
    noise_peak = abs(peaks(h));
    M = max(noise_peak(:));
    m = min(noise_peak(:));
    noise_peak = (noise_peak-m)./(M-m);
    Par.peak = noise_peak; 
    % add the spatially variant noise
%     Par.nim = zeros([h, w, ch]); % add noise
%     randn('seed',0);
%     for c = 1:ch
%         noise_this_channel = Sig_Level(c).*noise_peak.*randn([h, w]);
%         Par.nim(:, :, c) = Par.I(:, :, c) + noise_this_channel;
%     end
    % for a fair comparison, all the methods read the corrupted observations
    Par.nim = double( imread(fullfile(Noi_image_dir, Noi_im_dir(i).name)) );
    Par.nSig = Sig_Level .* mean( noise_peak(:) ); % 
    fprintf('%s c=%.2f delta=%.2f win=%d N=%d\nlambda=%.2f rho=%.2f alpha=%.2f t=%d\n', GT_im_dir(i).name(6:7), Par.Constant, Par.delta, Par.win, Par.nlsp, Par.lambda, Par.rho, Par.alpha, Par.t);
    PSNR0  =   csnr( Par.nim, Par.I, 0, 0 );
    SSIM0  =   cal_ssim( Par.nim, Par.I, 0, 0 );
    fprintf('  PSNR_0 = %f  SSIM_0 = %f\n', PSNR0, SSIM0);
    time0 = clock;
    [im_out, psnr_f, ssim_f] = TNF_Denoising_sv( Par.nim, Par.I, Par );
    PSNRs(i) = psnr_f; SSIMs(i) = ssim_f;
    Times(i) = etime(clock, time0);
    fprintf('Time used %f s\n', Times(i));
end
end
end
% fprintf('\n');
end
end