% script to save *.mrcs particle image stack to one *.mrcs stack file

path_particles = 'G:\20150205_sdp\Particles\SumCorr\';
pathout = 'G:\20150205_sdp\';

addpath(fullfile(cd, '../src/mrc'));

ds = 3;
vox_size = 1.32;

listing = dir([path_particles '*.mrcs']);
ntot = size(listing, 1);
ndow = round(ntot / ds);

for i = 1 : ndow
    fname = [path_particles '\' listing(i).name];
    a = ReadMRC(fname);
    if (i == 1)
        t = a;
    else
        t = cat(3, t, a);
    end
    
    perc = floor(i/ndow)*100;
    if (mod(i,round(0.2*ndow)) == 0)
        fprintf('%i', perc);
    elseif (mod(i,round(0.05*ndow)) == 0) 
        fprintf('.');
    end
end
fprintf('\n');
writeMRC(t, vox_size, [pathout '\' 'stackfile_ds' num2str(ds) '.mrcs']);
fprintf('done.\n');