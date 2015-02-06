function mask=get_mask_struct_ncd(xsize,pixtoedge)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Returns a mask structure of linear dim xsize(1)
%  Works only in 3D and 2D 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mask=struct('bin',[],'r',0);

if length(xsize)==3
    % 3D
    [x,y,z]=ndgrid(1:xsize(1));
    x=(x-(xsize(1)+1)/2);
    y=(y-(xsize(1)+1)/2);
    z=(z-(xsize(1)+1)/2);
    r=sqrt(x.^2+y.^2+z.^2);
elseif length(xsize)==2
    [x,y]=ndgrid(1:xsize(1));
    x=(x-(xsize(1)+1)/2);
    y=(y-(xsize(1)+1)/2);
    r=sqrt(x.^2+y.^2);
end

mask.r=(xsize(1)+1)/2-pixtoedge-0.5;
mask.bin= double(r<= mask.r);
    