function f_ind = get_freq_axis(voxel_size,im_size)

Fs = (1/voxel_size)/1e-10;
num = size(1:im_size/2,2);
f_ind = 1./(Fs/2*linspace(0,1,num));
f_ind = f_ind./1e-10;