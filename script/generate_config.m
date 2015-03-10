function outfile = generate_config(pathout, imfile, ctffile, substep, reconhalf, reconstartind, ...
    structfile, coordfile, maskfile, pf, pixfromedge, maxmem, numthreads, dispflag, ...
    f, numruns, maxnumiter, convtol, normprojint, rotstart, rotstep, rotend, transmax, transdelta, transwidth, ...
    alignims)
%GENERATE_CONFIG Generates txt config file given parameters
%   Run from cryo3d.m

fname = '_config_bm.txt';
fid = fopen([pathout fname], 'w');
fprintf(fid, 'imfile = %s \n', imfile);
if (~isempty(ctffile))
    fprintf(fid, 'ctffile = %s \n', ctffile);
end
fprintf(fid, 'substep = %i \n', substep);
fprintf(fid, 'reconhalf = %i \n', reconhalf);
fprintf(fid, 'reconstartind = %i \n', reconstartind);
fprintf(fid, 'structfile = %s \n', structfile);
fprintf(fid, 'coordfile = %s \n', coordfile);
if (~isempty(maskfile))
    fprintf(fid, 'maskfile = %s \n', maskfile);
end
fprintf(fid, 'pf = %i \n', pf);
fprintf(fid, 'pixfromedge = %i \n', pixfromedge);
fprintf(fid, 'maxmem = %i \n', maxmem);
fprintf(fid, 'numthreads = %i \n', numthreads);
fprintf(fid, 'dispflag = %i \n', dispflag);
fprintf(fid, 'f = %f \n', f);
fprintf(fid, 'numruns = %i \n', numruns);
fprintf(fid, 'maxnumiter = %i \n', maxnumiter);
fprintf(fid, 'convtol = %f \n', convtol);
fprintf(fid, 'normprojint = %i \n', normprojint);
fprintf(fid, 'rotstart = %f \n', rotstart);
fprintf(fid, 'rotstep = %f \n', rotstep);
fprintf(fid, 'rotend = %f \n', rotend);
fprintf(fid, 'transmax = %i \n', transmax);
fprintf(fid, 'transdelta = %i \n', transdelta);
fprintf(fid, 'transwidth = %i \n', transwidth);
fprintf(fid, 'alignims = %i \n', alignims);

fclose(fid);
outfile = [pathout fname];
end

