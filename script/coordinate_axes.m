% Function to create projection coordinates from given input
% Required parameters, example:
% pathout = 'G:\workspace\';
% dtheta = 5; Angular sampling in degrees
% Optional parameters, example:
% disp_flag = 1; Flag (0/1) for displaying the coordinates
% fig_num = 10; Figure number for displaying coordinates

function outfile = coordinate_axes(pathout, dtheta, disp_flag, fig_num)

if nargin < 4
    fig_num = 1;
end
if nargin < 3
    disp_flag = 0;
end
if nargin < 2
    disp('ERROR: Not enough input parameters');
    return;
end

addpath(fullfile(cd, '../src/best_match'));
coord_axes = create_approx_uniform_axes_whole_sphere(dtheta,disp_flag,fig_num);
outfile = [pathout 'coord_axes_' num2str(dtheta) '.mat'];
save(outfile, 'coord_axes');
