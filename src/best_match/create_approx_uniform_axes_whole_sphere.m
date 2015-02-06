function coord_axes=create_approx_uniform_axes_whole_sphere(dtheta,disp_flag,fig_num)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Creates coord axes with z-axis uniformly spaced on a polar grid
%   ntheta: number of zenith directions
%   nphi: number of azimuth direction
%   disp_flag: display flag, when 1 the normals are displayed in
%   fig(fig_num)
%   coord_axes: vector of projection axes
%   Always contains the vertical direction as a the first
%   normal.
%   Code history:
%       Modified from create_uiform_axes by H. Tagare  08/11/10.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

coord_axes=[1 0 0 0 1 0 0 0 1]';
alpha=pi/180;

for theta=dtheta:dtheta:180-dtheta,
    %theta
    dphi = dtheta / sind(theta);
    for phi=0:dphi:360-0.1,
        z=[sin(theta*alpha)*cos(phi*alpha) ...
            sin(theta*alpha)*sin(phi*alpha) ...
            cos(theta*alpha) ];
        tmp=get_axes_from_n_rot(z,0*alpha);
        coord_axes=[coord_axes tmp];
    end
end
last = [-1 0 0 0 1 0 0 0 -1]';
coord_axes = [coord_axes last];

if disp_flag==1
    figure(fig_num);
    plot_vectors(coord_axes);
end


function plot_vectors(coord_axes)
vsize=size(coord_axes);
for i=1:vsize(2),
        plot3(coord_axes(7,i),coord_axes(8,i),coord_axes(9,i),'k*');
        hold on;
        plot3([0 coord_axes(7,i)],[0 coord_axes(8,i)],[0 coord_axes(9,i)],'k');
end
plot3(0,0,0,'g+');
pause(0.1)