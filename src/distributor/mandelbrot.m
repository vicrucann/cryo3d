function count = mandelbrot( fname, resfname )
%MANDELBROT Mndelbrodt set calculation for given parameters
%   To use as a part of rshell-mat - bash script that helps to parallelize
%   matlab big data processing
%   2015 vicrucann@gmail.com

fprintf('The data file provided: %s\n', fname);

load(fname);
xGrid=xi;
yGrid=yi;
maxIterations=iter;
%------------------------

z0 = xGrid + 1i*yGrid;
count = ones( size(z0) );

z = z0;
for n = 0:maxIterations
    z = z.*z + z0;
    inside = abs( z )<=2;
    count = count + inside;
end
count = log( count );
save(resfname, 'count');
end

