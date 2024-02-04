function       CurPat = Image2Patch( im_out, par )
im_out     =  single(im_out);
CurPat     =  zeros(par.ps2ch, par.maxrc, 'double'); % 3p^2 * All
k          =  0;
for channel = 1:par.ch
    for i = 1:par.ps
        for j = 1:par.ps
            k    =  k+1;
            blk  =  im_out(i:end-par.ps+i,j:end-par.ps+j, channel);
            CurPat(k,:) = blk(:)';
        end
    end
end
