function mask=getVolMask(meanVol,thresh)
%Set up the filter
fsize=11;
[x y z]=meshgrid(1:fsize,1:fsize,1:fsize);
sigma=2;
c=fsize/2+0.5;
r2=((x-c).^2+(y-c).^2+(z-c).^2);
h=exp(-r2/(2*sigma^2));
h=h/sum(h(:));

%threshold the volume
mask=double(meanVol>=thresh);
%plot_surface(mask,0.1);

for k=1:2,
   mask=convn(mask,h,'same');
   mask=double(mask>=0.1);
    %plot_surface(mask,0.5);
end
mask=convn(mask,h,'same'); 


