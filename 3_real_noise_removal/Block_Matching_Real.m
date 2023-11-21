function       NL_mat = Block_Matching_Real( X, par)	% X := CurPat	[3p^2 * All]
% record the indexs of patches similar to the seed patch
NL_mat   =  zeros(par.nlsp, par.lenrc, 'single');

for  i  =  1 : par.lenrc  % 对于每一个key patch
    seed = X(:, par.SelfIndex(i));
        % seed[108,1] : 第i块key patch 的灰度值(3*36个位置)
    neighbor = X(:, par.NeighborIndex(1:par.NumIndex(i), i));
        % neighbor[108,4*sw^2] := 第i块key patch的4*sw^2个邻居的灰度值
    dis = sum(bsxfun(@minus, neighbor, seed).^2, 1);
        % neighbor的每一列(108,1)都减去seed再平方，再纵向求和
        % 即 计算每个邻居与这个key patch的相似度 (求和越大，越不相似)
        % dis 的第一列全0，因为这一列就是key patch本身
    [~,ind]   =  sort(dis);
        % 将dis升序排序，不破坏 dis
        % ind(i) := 排名第i的( 第i相似的patch )在dis中的下标
    indc = par.NeighborIndex( ind( 1:par.nlsp ), i );
        % 选取前 nlsp 个最相似的patch的下标, 组成(nlsp,1)的向量
    indc(indc == par.SelfIndex(i)) = indc(1); % added on 08/01/2017
    indc(1) = par.SelfIndex(i); % to make sure the first one of indc equals to off
    NL_mat(:, i) = indc;
end