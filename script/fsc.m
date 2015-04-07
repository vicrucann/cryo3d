% script to perform Fourier Shell Correlation

addpath(fullfile(cd, '../src/fsc'));
addpath(fullfile(cd, '../src/mrc'));
addpath(fullfile(cd, '../src/best_match'));

downsample = 2;
vol1 = ReadMRC('G:\workspace\db-hongwei\dselected_12\theta6\fbm_recon_masked_odd.mrc');
vol2 = ReadMRC('G:\workspace\db-hongwei\dselected_12\theta6\mask_resample_even.mrc');
vox_size = downsample*1.32;

fsc_with_mask(vol1, vol2, vox_size, 1, 6, 6);