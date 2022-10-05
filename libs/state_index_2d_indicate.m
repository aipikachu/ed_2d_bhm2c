function [idx_st] = state_index_2d_indicate(st_up,st_dn,basis)


%%
nMax = basis.nMax;
L = basis.L;
idxlt = basis.idxstatel;


%%
stat_up_to1d_chain = reshape(st_up',1,[]);
stat_dn_to1d_chain = reshape(st_dn',1,[]);

idx_up_cur = stat_up_to1d_chain * ((nMax+1).^(L-1:-1:0))';
idx_dn_cur = stat_dn_to1d_chain * ((nMax+1).^(L-1:-1:0))';


%%
[~,idx_st] = ismember([idx_up_cur,idx_dn_cur],idxlt,'rows');

