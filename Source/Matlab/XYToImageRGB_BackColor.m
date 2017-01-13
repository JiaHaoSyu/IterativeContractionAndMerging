function [I_ori,I_skewed,Scatter_map] = XYToImageRGB_BackColor(data,XY,BackColor)
% syntax:
%       [I_ori,I_skewed,Scatter_map] = XYToImageRGB_BackColor(data,XY,BackColor)
% description:
%       Display contracted image such as Fig.4(b)-(d) with specific background
%
% Input arguments:
% data: original input image with MxNx3 matrix
% XY : Twisted spatial coordinate with MxNx2 matrix         
% BackColor: background color with 1x3 matrix

% Output arguments:
% I_ori :  original input image
% I_skewed : contracted image
% Scatter_map : mask of Scatter_map
%  Jia-Hao Syu <jiahaushiu4@gmail.com>
%  December 2016

[img_h_ori, img_w_ori, img_c_ori] = size(data);
I_i = min(max(floor(XY(:,:,1)+0.5),1),img_h_ori);
I_j = min(max(floor(XY(:,:,2)+0.5),1),img_w_ori);
I_ind = sub2ind([img_h_ori, img_w_ori], I_i, I_j);
I_skewed = zeros(size(data));
Scatter_map = zeros([img_h_ori, img_w_ori]);
I_ori = reshape(data,[img_h_ori*img_w_ori, img_c_ori]);
I_skewed = reshape(I_skewed,[img_h_ori*img_w_ori, img_c_ori]);

I_indS = full(sparse(I_ind(:),ones(size(I_ind(:),1),1),ones(size(I_ind(:),1),1),size(I_ind(:),1),1));
Scatter_map = reshape(I_indS,[img_h_ori, img_w_ori]);

I_indR = full(sparse(I_ind(:),ones(size(I_ind(:),1),1),I_ori(:,1),size(I_ind(:),1),1));
I_indG = full(sparse(I_ind(:),ones(size(I_ind(:),1),1),I_ori(:,2),size(I_ind(:),1),1));
I_indB = full(sparse(I_ind(:),ones(size(I_ind(:),1),1),I_ori(:,3),size(I_ind(:),1),1));

I_skewed(:,1) = I_indR;I_skewed(:,2) = I_indG;I_skewed(:,3) = I_indB;

for k = 1:img_c_ori
  I_skewed(Scatter_map>0,k) = I_skewed(Scatter_map>0,k)./Scatter_map(Scatter_map>0);
  I_skewed(Scatter_map==0,k) = BackColor(k);
end

I_skewed = reshape(I_skewed,[img_h_ori,img_w_ori,img_c_ori]);