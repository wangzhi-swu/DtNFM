clear; 

Original_image_dir  =    'kodak_color/'; 
fpath = fullfile(Original_image_dir, '*.png');
im_dir  = dir(fpath);
im_num = length(im_dir);

nSig = [20 35 5];
% nSig = [30 10 50];

PSNR_arr = zeros([im_num,1]);
SSIM_arr = zeros([im_num,1]);
Time_arr = zeros([im_num,1]);
for ALPHA = [1.8 ]
for Truncate=[2]
for LAMBDA = [0.8]
for RHO = [0.3]
for i = 1:24
Par.nSig      =   nSig;       % STD of the noise image
Par.win       =   20;         % Non-local patch searching window
Par.Constant  =   8*sqrt(2);  % Constant num for the weight vector 
Par.Innerloop =   2;          % InnerLoop Num of between re-blockmatching
Par.ps        =   6;          % Patch size, larger values would get better performance, but will be slower
Par.step      =   5;          % The step between neighbor patches, smaller values would get better performance, but will be slower
Par.Iter      =   12; 
Par.display = true;
Par.delta     =   0.1; % step back length of iterative regularization
Par.maxIter   =   10;

Par.rho       =   RHO;    
Par.mu        =   1.002;

Par.lambda    =   LAMBDA;
Par.t = Truncate; 

    Par.alpha = ALPHA;
    Par.image   =   i;
    Par.nSig0   =   nSig;
    Par.nlsp    =   70; 
    Par.I = double( imread(fullfile(Original_image_dir, im_dir(i).name)) );
    S = regexp(im_dir(i).name, '\.', 'split');
    [h, w, ch] = size(Par.I);

    Par.nim = zeros([h, w, ch]); % add noise
    randn('seed',0);
    for c = 1:ch
        Par.nim(:, :, c) = Par.I(:, :, c) + Par.nSig0(c) * randn([h, w]);
    end

    fprintf('%s c=%.2f lambda=%.2f rho=%.2f alpha=%.2f t=%d N=%d\n', im_dir(i).name(6:7), Par.Constant, Par.lambda, Par.rho, Par.alpha, Par.t, Par.nlsp);
    PSNR0  =   csnr( Par.nim, Par.I, 0, 0 );
    SSIM0  =   cal_ssim( Par.nim, Par.I, 0, 0 );

    fprintf('PSNR_0 = %f  SSIM_0 = %f\n', PSNR0, SSIM0);

    time0 = clock;
    [im_out, Par, psnr_final, ssim_final] = TNF_Denoising( Par.nim, Par.I, Par );
    Time_arr(i) = etime(clock, time0);
    fprintf('Time used %f s\n', Time_arr(i));
    PSNR_arr(i) = psnr_final; SSIM_arr(i) = ssim_final;
end
end
end
end
end
