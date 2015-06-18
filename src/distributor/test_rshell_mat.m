%% Example of matlab scrips that launches bash script to control data parallelization among servers
% One examples is considered: calculation of madelbrot set

%% Setting up
clc; clear; close all;
login = 'cryo';
ppath = '/home/cryo/dop'; % distributed operations, destination on remote
ipaddrs = ['130.132.104.236' ' ' '172.21.9.92' ' ' '172.23.2.105' ' ' '172.23.5.77']; % list of ip addresses
remmat = 'mandelbrot'; % name of matlab function that will be launched on remote server
varmat = 'mnd'; % when splitting data, they will be saved under varmat.mat name on disk
sleeptime = 5;
resfold = 'dres'; % name of the result folder
bashscript = fullfile(pwd,'dhead.sh'); % main bash script that organizes data processing
printout = 0; % print the bash output (1) or not (0)

[ncluster ~] = find(ipaddrs==' '); % to break data into n clusters (as many as given servers)
ncluster = size(ncluster,2)+1;

%% Mandelbrot set
% Given resolution and iteration number, find corresponding Mandelbrot set

% input parameters
iter = 2000;
isize = 4000;
xlim = [-2, 1]; % to split
ylim = [-1.5, 1.5]; % to split
figpos = [100 100 1000 1000];

% pre-calcualtion
x = linspace( xlim(1), xlim(2), isize );
y = linspace( ylim(1), ylim(2), isize );
[xGrid,yGrid] = meshgrid( x, y );
szx = ceil(size(xGrid,2)/ncluster);

% perform calculation by distributing among the servers
% split (assume we break along "X" dimension)
fprintf('\n\nCalculation using remotes \n')
tic();
for i=1:ncluster
    if (i ~= ncluster)
        xi=xGrid(:, szx*(i-1)+1:szx*i);
        yi=yGrid(:, szx*(i-1)+1:szx*i);
    else
        xi=xGrid(:, szx*(i-1)+1:end);
        yi=yGrid(:, szx*(i-1)+1:end);
    end
    save([varmat int2str(i) '.mat'], 'xi', 'yi', 'iter');
end
system(['chmod u+x ' bashscript])
if printout
    cmdStr = [bashscript ' ' login ' ' ppath ' ' ipaddrs ' ' remmat ' ' varmat ' ' int2str(sleeptime) ' ' resfold];
else
    cmdStr = [bashscript ' ' login ' ' ppath ' ' ipaddrs ' '...
        remmat ' ' varmat ' ' int2str(sleeptime) ' ' resfold '>' remmat '.log 2>&1'];
end
% perform the command
system(cmdStr)

% merge the results
res = zeros(isize, isize);
for i=1:ncluster
    load([resfold '/' 'result_' varmat int2str(i) '.mat']);
    if (i ~= ncluster)
        res(:,szx*(i-1)+1:szx*i) = count;
    else
        res(:,szx*(i-1)+1:end) = count;
    end
end
tcluster = toc();

% perform the full calculation of mandelbrot on local
fprintf('Calculation on local...\n');
tic();
z0 = xGrid + 1i*yGrid;
count0 = ones( size(z0) );
z = z0;
for n = 0:iter
    z = z.*z + z0;
    inside = abs( z )<=2;
    count0 = count0 + inside;
    if (mod(n, iter*0.25) == 0)
        fprintf('%i', ceil(n/iter*100));
    elseif (mod(n,iter*0.05) == 0)
        fprintf('.');
    end
end
count0 = log( count0 );
tlocal=toc();
fprintf( ' -> done\n');

fprintf('\n\nLocal time vs distributed time: \n    %1.2fsecs vs %1.2fsecs \n', tlocal, tcluster);

% display local result
figure;
fig = gcf;
fig.Position = figpos;
imagesc( x, y, count0 );
axis image
colormap( [jet();flipud( jet() );0 0 0] );

% display distributed result
figure;
fig = gcf;
fig.Position = figpos;
imagesc( x, y, res );
axis image
colormap( [jet();flipud( jet() );0 0 0] );

