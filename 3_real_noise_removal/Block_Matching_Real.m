function       NL_mat = Block_Matching_Real( X, par)	% X := CurPat	[3p^2 * All]
% record the indexs of patches similar to the seed patch
NL_mat   =  zeros(par.nlsp, par.lenrc, 'single');

for  i  =  1 : par.lenrc  % ����ÿһ��key patch
    seed = X(:, par.SelfIndex(i));
        % seed[108,1] : ��i��key patch �ĻҶ�ֵ(3*36��λ��)
    neighbor = X(:, par.NeighborIndex(1:par.NumIndex(i), i));
        % neighbor[108,4*sw^2] := ��i��key patch��4*sw^2���ھӵĻҶ�ֵ
    dis = sum(bsxfun(@minus, neighbor, seed).^2, 1);
        % neighbor��ÿһ��(108,1)����ȥseed��ƽ�������������
        % �� ����ÿ���ھ������key patch�����ƶ� (���Խ��Խ������)
        % dis �ĵ�һ��ȫ0����Ϊ��һ�о���key patch����
    [~,ind]   =  sort(dis);
        % ��dis�������򣬲��ƻ� dis
        % ind(i) := ������i��( ��i���Ƶ�patch )��dis�е��±�
    indc = par.NeighborIndex( ind( 1:par.nlsp ), i );
        % ѡȡǰ nlsp �������Ƶ�patch���±�, ���(nlsp,1)������
    indc(indc == par.SelfIndex(i)) = indc(1); % added on 08/01/2017
    indc(1) = par.SelfIndex(i); % to make sure the first one of indc equals to off
    NL_mat(:, i) = indc;
end