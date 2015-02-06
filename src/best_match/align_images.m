function [ aligned_ims ] = align_images(ims,rotinds,transinds,rots,trans)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

numim = size(ims,3);
n = size(ims,1);
aligned_ims = zeros(size(ims),'single');
for i = 1:numim
    dx = -trans(transinds(i),1);
    dy = -trans(transinds(i),2);
    transim = zeros(n,n,'single');
    if dy < 0
        if dx < 0
            transim(1:end+dy,1:end+dx) = ims(1-dy:end,1-dx:end,i);
        else
            transim(1:end+dy,1+dx:end) = ims(1-dy:end,1:end-dx,i);
        end
    else
        if dx < 0
            transim(1+dy:end,1:end+dx) = ims(1:end-dy,1-dx:end,i);
        else
            transim(1+dy:end,1+dx:end) = ims(1:end-dy,1:end-dx,i);
        end
    end
    aligned_ims(:,:,i) = imrotate(transim,-rots(rotinds(i)),'bilinear','crop');
end


end

