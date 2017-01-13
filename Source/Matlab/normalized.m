function [Input] = normalized(Input)
[~,~,dim] = size(Input);

for i = 1:dim
    Input(:,:,i) = (Input(:,:,i)-min(min(Input(:,:,i))))/(max(max(Input(:,:,i)))-min(min(Input(:,:,i))));
end