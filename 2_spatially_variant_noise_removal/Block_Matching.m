function       NL_mat = Block_Matching( X, par) % X ~ CurPat[108, 386841]，blk_arr ~ NL_mat
% record the indexs of patches similar to the seed patch
NL_mat   =  zeros(par.nlsp, par.lenrc, 'single');

for  i  =  1 : par.lenrc  % for each key patch
    seed = X(:, par.SelfIndex(i));
    neighbor = X(:, par.NeighborIndex(1:par.NumIndex(i), i));
    dis = sum(bsxfun(@minus, neighbor, seed).^2, 1);
    [~,ind]   =  sort(dis);
    indc = par.NeighborIndex( ind( 1:par.nlsp ), i );
    NL_mat(:, i) = indc;
end