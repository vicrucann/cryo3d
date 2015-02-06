function T=get_axes_from_n_rot(n,phi_r)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Returns the xyz axes using normal 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


proj_dir=n;
%proj_dir=proj_dir/sqrt(proj_dir*proj_dir');  %the projection direction
z=[0 0 1]; %The z axis
%Check if the proj direction is z
a=cross(z,proj_dir);
if (sqrt(sum(a.*a)) < 0.0001)
    %Projection direction is z
    T1=eye(3,3);
    T2=eye(3,3);
else
    %Projection direction is not z
    %Rotate x-y plane
    phi=atan2(proj_dir(2),proj_dir(1));
    T1=[cos(phi) -sin(phi) 0;...
        sin(phi) cos(phi) 0;...
        0 0 1];

    %Rotate around y axis
    theta=atan2(sqrt(proj_dir(1)^2+proj_dir(2)^2),proj_dir(3));
    T2=[cos(theta) 0 sin(theta);...
        0 1 0;...
        -sin(theta) 0 cos(theta)];
end
%Final in plane xy rotation
T3=[cos(phi_r) -sin(phi_r) 0;...
        sin(phi_r) cos(phi_r) 0;...
        0 0 1];
    
%Net rotation
T=(T1*T2)*T3;
T=[T(:,1); T(:,2); T(:,3)];

