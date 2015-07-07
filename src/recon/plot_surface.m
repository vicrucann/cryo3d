function plot_surface(x,t)
clf;
p = patch(isosurface(x,t));
axis([0 size(x,1) 0 size(x,2) 0 size(x,3)]);
set(p, 'FaceColor', 'red', 'EdgeColor', 'none');
daspect([1 1 1]);
view(3)
camlight
lighting gouraud
drawnow;
