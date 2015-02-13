function passed = coordinate_axes(dtheta,disp_flag,fig_num)

passed = 0;
addpath(fullfile(cd, '../src/best_match'));
coord_axes = create_approx_uniform_axes_whole_sphere(dtheta,disp_flag,fig_num);

passed = 1;