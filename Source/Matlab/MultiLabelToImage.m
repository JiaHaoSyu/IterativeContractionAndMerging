function [LabelImage] = MultiLabelToImage(LabelCombine)
% syntax:
%       [LabelImage] = MultiLabelToImage(LabelCombine)
% description:
%       Multiple Labeling image(segments) to Multiple color image with
%       consistent color change
%
% Input arguments:
    % LabelCombine : Multiple Labeling image(segments) with MxNxK matrix, K
    % indicates the K multiple Labeling images(segments)
   
% Output arguments:
    % LabelImage : Multiple color image with MxNx3xK matrix

%   Jia-Hao Syu <jiahaushiu4@gmail.com>
%  December 2016

[~,~,dim] = size(LabelCombine);

for i = dim:-1:1
    if i==dim % 1 segments
        Label = LabelCombine(:,:,i);
        R = sparse(ones(max(Label(:)),1),1:max(Label(:)),rand([max(Label(:)) 1]));
        G = sparse(ones(max(Label(:)),1),1:max(Label(:)),rand([max(Label(:)) 1]));
        B = sparse(ones(max(Label(:)),1),1:max(Label(:)),rand([max(Label(:)) 1]));
              
        OutputImg = full(R(Label));
        OutputImg(:,:,2) = full(G(Label));
        OutputImg(:,:,3) = full(B(Label));
        LabelImage(:,:,:,i) = OutputImg;
    else
        LabelTarge = Label;
        LabelNew = LabelCombine(:,:,i);
        % Larger overlap area provides same color, others give random color
        Map = sparse(LabelTarge(:),LabelNew(:),ones(size(Label(:))));
        [~,ia] = max(Map,[],2);
        
        [cc,ic,~] = unique(ia,'sorted'); % cc is unique labeling
%         [cc1 ic1 ci1] = unique(id);
        % non-Overlap area : give a random color
        
        if verLessThan('matlab', '8.0')
            [bb,~,~] = setxor(cc,[1:max(LabelNew(:))]);
            % -- Put code to run under MATLAB 8.0 and earlier here --
        else
            [bb,~,~] = setxor(cc,[1:max(LabelNew(:))],'legacy');
            % -- Put code to run under MATLAB 8.01 and later here --
        end
                               
        Ra = rand([size(bb,2) 1]);
        Ga = rand([size(bb,2) 1]);
        Ba = rand([size(bb,2) 1]);
        
        R1 = full(R(ic));
        G1 = full(G(ic));
        B1 = full(B(ic));
        RNew = [R1(:);Ra];
        GNew = [G1(:);Ga]; 
        BNew = [B1(:);Ba]; 
        row_inds = [cc;bb'];
        col_inds = ones(max(LabelNew(:)),1);
             
        R = sparse(row_inds(:),col_inds(:),RNew(:));
        G = sparse(row_inds(:),col_inds(:),GNew(:));
        B = sparse(row_inds(:),col_inds(:),BNew(:));
        OutputImg = full(R(LabelNew));
        OutputImg(:,:,2) = full(G(LabelNew));
        OutputImg(:,:,3) = full(B(LabelNew));
        LabelImage(:,:,:,i) = OutputImg;
        Label = LabelNew;
    end
end


