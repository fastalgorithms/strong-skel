addpath('../fortran_src')
addpath('../')
addpath('../sv')
addpath('../mv')
addpath('../src')
addpath('../src/helm_dirichlet/');
addpath('../src/helm_sound_hard/');
run('../../FLAM/startup.m');

% Patch occupancy
occ = 50;

% Number of proxy points, not used?
nproxy = 50;

% Final tolerance for factorisation
rank_or_tol = 0.51e-4;

% Number of patches in u dir on torus
npu = 10;

% Related to number of quadrature points
norder = 3;
zk = 1.0;

% Initialize a wiggly torus
radii = [1;3;0.25]; 
scales = [1.2;1.0;1.7];

nnu = npu;
nnv = npu;
nosc = 5;
sinfo = wtorus(radii, scales, nosc, nnu, nnv, norder);

x = sinfo.srcvals(1:3,:);

shp = alphaShape(x(1, :)', x(2,:)', x(3,:)', 1.5,"HoleThreshold",0);
plot(shp) 
 
[tri, xyz] = boundaryFacets(shp);

h = trisurf(tri, xyz(:,1), xyz(:,2), xyz(:,3), 'FaceColor', 'interp', 'FaceAlpha',1);

h.LineStyle='none';

colormap('copper');

l = light('Position',[-1 0 0], 'Style','infinite');
lighting gouraud;
material shiny;
