clear; 

Original_image_dir  =    'kodak_color/'; 
fpath = fullfile(Original_image_dir, '*.png');
im_dir  = dir(fpath);
im_num = length(im_dir);

% nSig = [20 35 5];
nSig = [30 10 50];

for ALPHA = [3]
for Truncate=[3]
for LAMBDA = [1]
% lb = -LAMBDA + 1.5;
PSNRs = zeros([im_num, 1]); SSIMs = zeros([im_num, 1]); Times = zeros([im_num, 1]);cnt=0;
for RHO = [0.6]%[lb : 0.05: lb + 0.11]
cnt=cnt+1;
Par.nSig      =   nSig;          % STD of the noise image
Par.win       =   20;            % Non-local patch searching window
Par.Constant  =   6*sqrt(2);   % Constant num for the weight vector
Par.Innerloop =   2;             % InnerLoop Num of between re-blockmatching
Par.ps        =   6;             % Patch size, larger values would get better performance, but will be slower
Par.step      =   5;             % The step between neighbor patches, smaller values would get better performance, but will be slower
Par.Iter      =   12;            % total iter numbers
Par.display = true;
Par.delta     =   0.1;           % iterative regularization parameter
Par.maxIter   =   10;

Par.rho       =   RHO;           % In final version, this parameter will be changed
Par.mu        =   1.002;

% this parameter is not finally determined yet
Par.lambda    =   LAMBDA;        % for different noise levels, this parameter should be tuned to achieve better performance
Par.t = Truncate; 
for i = 19
% record all the results in each iteration
Par.PSNR = zeros(Par.Iter, im_num, 'single');
Par.SSIM = zeros(Par.Iter, im_num, 'single');
    Par.alpha = ALPHA;
    Par.image   =   i;
    Par.nSig0   =   nSig;
    Par.nlsp    =   80;   % Initial Non-local Similar Patches number
    Par.I = double( imread(fullfile(Original_image_dir, im_dir(i).name)) );
    S = regexp(im_dir(i).name, '\.', 'split');
    [h, w, ch] = size(Par.I);

    Par.nim = zeros([h, w, ch]); % add noise
    randn('seed',0);
    for c = 1:ch
        Par.nim(:, :, c) = Par.I(:, :, c) + Par.nSig0(c) * randn([h, w]);
    end

    fprintf('%s c=%.2f lambda=%.2f rho=%.2f alpha=%.2f t=%d N=%d\n', im_dir(i).name(6:7), Par.Constant, Par.lambda, Par.rho, Par.alpha, Par.t, Par.nlsp);
    PSNR  =   csnr( Par.nim, Par.I, 0, 0 );
    SSIM  =   cal_ssim( Par.nim, Par.I, 0, 0 );

    fprintf('PSNR_0 = %f  SSIM_0 = %f\n', PSNR, SSIM);

    time0 = clock;
    [im_out, psnr_f, ssim_f] = TNF_Denoising( Par.nim, Par.I, Par );
    Times(i,cnt) = etime(clock, time0);
    fprintf('Time used %f s\n',Times(i,cnt));
    PSNRs(i,cnt) = psnr_f; SSIMs(i,cnt) = ssim_f;
end

end
end
% fprintf('\n\n');
end
end