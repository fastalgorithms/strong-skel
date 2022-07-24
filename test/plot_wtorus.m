addpath('../fortran_src')
addpath('../')
addpath('../sv')
addpath('../mv')
addpath('../src')
addpath('../src/helm_dirichlet/');
addpath('../src/helm_sound_hard/');
run('../../FLAM/startup.m');

uu = linspace(0, 2*pi, 100);
vv = linspace(0, 2*pi, 50);
rr = 1;
[Uu, Vv] = meshgrid(uu, vv);

x = (rr.*cos(Uu) + 2 + 0.25*cos(5*Vv)).*cos(Vv)*1.2;
y = (rr.*cos(Uu) + 2 + 0.25*cos(5*Vv)).*sin(Vv)*1.0;
z = rr.*sin(Uu)*1.7;

h = surf(x,y,z);
h.LineStyle='none';

colormap('copper');

l = light('Position',[-1 0 0], 'Style','infinite');
lighting gouraud;
material shiny;

xlabel('x')
% set(gca, "Visible", 'off')
% exportgraphics(h, 'wiggly_torus.eps', 'Resolution', '300', 'ContentType', 'vector')