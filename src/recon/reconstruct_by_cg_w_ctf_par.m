function [x_r, dx, iter] =reconstruct_by_cg_w_ctf_par(proj,data_axes,ctfs,mask,l_norm,l_smooth,iter_lim,stop_lim,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Reconstructs a volume by spatial conjugte gradient
%   proj: N X N X K array of N X N images along K projection directions
%   normals: K X 3 array of projection directions
%   ctf_params: M X 2 array of ctf parameters
%   mask: N X N X N 0-1 array defining the mask within which x is
%   reconstructed
%   l_norm: Norm regularization constant
%   l_smooth: Smoothness regularization constant
%   iter_lim: iteration limit (max no. of iterations)
%   stop_lim: fractional change in reconst taken as stopping criterion
%   x_r: reconstructed volume N X N X N
%   This function needs CG_UTIL environmental variable to be set
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Original code by Hemant Tagare
% Modified by Nicha C. Dvornek, 11/2013 for parallel processing

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Set the utilities directory
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Get the utilities directory
% util_dir=getenv('CG_UTIL');
% if isempty(util_dir)
%     disp('reconstruct_by_cg: CG_UTIL environment variable not set');
%     return
% else
%     path(util_dir,path);
% end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Computations  begin here
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55

if length(varargin) < 1
    %Initialize the recontruction
    xsize=size(mask.bin);
    x_r=zeros(xsize);
else
    x_r = varargin{1};
end

%initialize algorithm parameters
eps_denom=1e-20;
stop_flag=0;
iter=0;
x_prev=x_r;
dx = inf;

%Outer iteration
while stop_flag==0
    %Get gradient and display norm

    g=get_grad_with_ctf_par(x_r,mask.r,proj,data_axes,l_norm,l_smooth,ctfs);


    %Get conjugate direction
    if iter==0
        %First iteration
         d= -g;
    else
        %Calculate beta
        beta=g(:)'*g(:)/(g_prev(:)'*g_prev(:)+eps_denom);
        %Update conjugate direction
        d= -g +beta*d_prev;
    end
    
    %Get step size
    alpha= -g(:)'*d(:)/( get_dQd_with_ctf_par(d,mask.r,data_axes,l_norm,l_smooth,ctfs)+eps_denom);
    
    %Update reconstruction
    x_r=(x_r+alpha*d).*mask.bin;
    
    %Check for stopping
    iter=iter+1;

    if iter >= iter_lim
        stop_flag=1;
    end
    dx_last = dx;
    x_diff=x_r-x_prev;
    x_diff_norm=sqrt(x_diff(:)'*x_diff(:));
    x_norm=sqrt(x_r(:)'*x_r(:));
    dx=x_diff_norm/(x_norm+eps_denom);
    if dx < stop_lim
        stop_flag=1;
    end
    if dx_last < dx && iter > 3
        stop_flag = 1;
        x_r = x_prev;
    end
    
    disp(sprintf('reconstruct_by_cg_w_ctf: iter = %d dx=%f\n',iter,dx));
    
    % Update prevs
    g_prev=g;
    d_prev=d;
    x_prev=x_r;
end
