function [Label,A1,A,C] = GraphAndConnectedComp(MapCombineAll,row_inds,col_inds,tvals)
% Merging with graph connected component
[I_Height,I_Width,~] = size(MapCombineAll);
A=sparse(row_inds,col_inds,tvals,I_Height*I_Width,I_Height*I_Width);
sumA=sum(A,2);
A1 = spdiags(sumA(:),0,I_Height*I_Width,I_Height*I_Width)-A;
[~,C] = conncomp(A);
Label = reshape(C,[I_Height,I_Width]);