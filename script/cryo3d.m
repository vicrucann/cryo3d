% script to run the whole pipeline
% all the parameters are set here

%% Most frequent edited parameters

pathout = 'G:\workspace\db-hongwei\ds03rd_lpf30_ds2\test\';
pathdata = 'G:\20150205_sdp\';
caching = 0;
dtheta = 12;
maskfile = [];
maxmem = 26624;
numthreads = 6;
rotstep = 5;
transmax = 4;

%% Preprocessing

structfile = [pathdata '\' 'run1_class001.mrc'];
stackfile = [pathdata '\' 'stackfile_ds3_third.mrcs'];
ctffile = [];
lpf = 30;
sigma = 1;
ds = 2; % downsampling of initial
downsample = 2;

structfile = init_volume(pathout, structfile, lpf, sigma, ds);
coordfile = coordinate_axes(pathout, dtheta);
imfile = preprocess_images(pathout, stackfile, ctffile, donwsample);

%% Best Match

configfile = '';
% Parameters for only using part of the data
substep = 0;
reconhalf = 0;
reconstartind = 0;
% Initial Model and Related Parameters
pf = 1;
pixfromedge = 4;
% System Parameters
dispflag = 0;
% Algorithm parameters
f = 0.0001;
numruns = 2;
maxnumiter = 10;
convtol = 0.01;
normprojint = 1;
% Transformation Search Parameters
rotstart = 1;
rotend = 359.9;
transdelta = 1;
transwidth = 1;
% Output parameters
alignims = 1;
% Generate config file

best_match(pathout, configfile, caching);