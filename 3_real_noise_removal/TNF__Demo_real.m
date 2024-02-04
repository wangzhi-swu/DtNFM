clear;
GT_Original_image_dir = 'CC15/';
GT_fpath   = fullfile(GT_Original_image_dir, '*mean.png');
TT_Original_image_dir = 'CC15/';
TT_fpath   = fullfile(TT_Original_image_dir, '*real.png');
GT_im_dir  = dir(GT_fpath);
TT_im_dir  = dir(TT_fpath);
im_num     = length(TT_im_dir);

Par.win       =   20; 
Par.delta     =   0; 
Par.Constant  =   16*sqrt(2);
Par.Innerloop =   2;
Par.ps        =   6; 
Par.step      =   5;
Par.Iter      =   12; 
Par.maxIter   =   10;

for i = [1]
for ALPHA = [2]
    Par.t = 0;                                                                                                                                                                                            
    Par.alpha = ALPHA;
for lambda = [2.3]
    Par.lambda = lambda;
    for rho = [0.9]
        Par.rho = rho;
        Par.PSNR = zeros(Par.Iter, im_num, 'single');
        Par.SSIM = zeros(Par.Iter, im_num, 'single');

        Par.mu = 1.002;
        Par.image = i;
        Par.nlsp  =  70;
        Par.I = double( imread(fullfile(GT_Original_image_dir, GT_im_dir(i).name)) );
        S = regexp(GT_im_dir(i).name, '\.', 'split');
        [h, w, ch] = size(Par.I);
        Par.nim = double( imread(fullfile(TT_Original_image_dir, TT_im_dir(i).name)) );
        for c = 1:ch
            Par.nSig0(c,1) = NoiseEstimation(Par.nim(:, :, c), Par.ps);
        end
        fprintf('#%d c=%.2f lmd=%.2f rho=%.2f alp=%.2f t=%d\n', i, Par.Constant, Par.lambda, Par.rho, Par.alpha, Par.t);
        PSNR = csnr( Par.nim, Par.I, 0, 0 );
        SSIM = cal_ssim( Par.nim, Par.I, 0, 0 );
        fprintf('  PSNR_0 = %f  SSIM_0 = %f\n', PSNR, SSIM);

        time0 = clock;
        [im_out, Par] = TNF_Denoising_real( Par.nim, Par.I, Par );
        fprintf('Total elapsed time = %f s\n', (etime(clock,time0)) );
        im_out(im_out>255) = 255;
        im_out(im_out<0.0) = 0.0;
    end
end
end
end