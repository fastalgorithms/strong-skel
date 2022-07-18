% function [varargout] =  ie_fmm3dbie_sound_hard(igeomtype,iref,npu,norder,ndir,zk)
% IE_SPHERE  An example usage of strong skeletonization, solving a
%  second-kind integral equation (Helmholtz combined field potential) on a
%  wiggly torus
%
%  - OCC:         The occupancy parameter, specifying the maximum number of 
%                 points a node in the octree can contain before it is
%                 subdivided.  This therefore gives an upper bound on the 
%                 number of points in a leaf node.
%  - P:           The number of proxy points to use to discretize the proxy 
%                 sphere used during skeletonization.
%  - RANK_OR_TOL: If a natural number, the maximum number of skeletons to
%                 select during a single box of skeletonization.  If a
%                 float between 0 and 1, an approximate relative tolerance
%          \       used to automatically select the number of skeletons.

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
npu = 5;

% Related to number of quadrature points?
norder = 3;
zk = 1.0;

% Initialize a wiggly torus
radii = [1.0;2.0;0.25]; 
scales = [1.2;1.0;1.7];

nnu = npu;
nnv = npu;
nosc = 5;
sinfo = wtorus(radii, scales, nosc, nnu, nnv, norder);

% Test points
m = 1;

% Seed for random numbers
rng(42);

% Pick 'm' points on object surface
xyz_in = zeros(3, m);
uu = rand(m,1)*2*pi;
vv = rand(m,1)*2*pi;
rr = rand(m,1)*0.67;

uu = 0.3;
vv = 0.5;
rr = 0.6;

uu = 0;
vv = 0;
rr = 0.01;

xyz_in(1,:) = (rr.*cos(uu) + 2 + 0.25*cos(5*vv)).*cos(vv)*1.2;
xyz_in(2,:) = (rr.*cos(uu) + 2 + 0.25*cos(5*vv)).*sin(vv)*1.0;
xyz_in(3,:) = rr.*sin(uu)*1.7;

% Pick 'm' points on parallel surface to the object
xyz_out = zeros(3, m);
uu = rand(m, 1)*2*pi;
vv = rand(m, 1)*2*pi;
rr = rand(m, 1)*0.67 + 1.33;



xyz_out(1,:) = (rr.*cos(uu) + 2 + 0.25*cos(5*vv)).*cos(vv)*1.2;
xyz_out(2,:) = (rr.*cos(uu) + 2 + 0.25*cos(5*vv)).*sin(vv)*1.0;
xyz_out(3,:) = rr.*sin(uu)*1.7;

xyz_out = [-3.5; 10; 5];

% Quadrature points
x = sinfo.srcvals(1:3,:);
x_or = x;
x = repmat(x, [1,2]);
nu = sinfo.srcvals(10:12,:);
area = sinfo.wts';
N = sinfo.npts;

% Complex parameterizations
zpars = complex([zk; -1j*zk; 1.0]);
zstmp = complex([zk;1.0;0.0]);
zdtmp = complex([zk;0.0;1.0]);
eps = rank_or_tol;

% Compute the quadrature corrections
tic, Q = helm_neu_near_corr(sinfo, zpars(1:2), eps); tquad=  toc;
S = Q{1};
P = zeros(N,1);
w = whos('S');
fprintf('quad: %10.4e (s) / %6.2f (MB)\n',tquad,w.bytes/1e6)

% Set system matrix
Afun_use = @(i, j) Afun_helm_sound_hard_wrapper(i, j, x_or, zpars, nu, area, P, S);

% Set proxy function
pxyfun_use = @(x, slf, nbr, proxy, l, ctr) pxyfun_helm_sound_hard_wrapper(x, slf, nbr, proxy, l, ctr, zpars, nu, area);

% Factorize the matrix
opts = struct('verb', 1, 'symm','n', 'zk', zk);
tic, F = srskelf_asym_new(Afun_use, x, occ, rank_or_tol, pxyfun_use, opts); tfac = toc;
% w = whos('F');
% fprintf([repmat('-',1,80) '\n'])
% fprintf('mem: %6.4f (GB)\n',w.bytes/1048576/1024)
% save('saveF.mat', 'F');
% Compute Neumann data - g = \Del_x S(x,y)
% 
% q = 1;
% x_or = sinfo.srcvals(1:3,:);
% B = -1*helm_sound_hard_kernel(x_or, xyz_in, zk, nu)*q;
% B = B.*sqrt(area).';
% 
% return

% A = Afun_use(1:2*N, 1:2*N);

% Load factorization from disk
% F = matfile('saveF.mat').F;

% Random weights
q = rand(m, 1)-0.5 + 1j*(rand(m,1)-0.5);
q = 1;

% Doesn't matter, as calculating slp on surface?
nu2 = zeros(3, m);

% Compute Neumann data - g = \Del_x S(x,y)
x_or = sinfo.srcvals(1:3,:);
B = -1*helm_sound_hard_kernel(x_or, xyz_in, zk, nu)*q;
B = B.*sqrt(area).';

% Add zero data due to system - [-g, 0, ..., -g, 0]
Btmp = zeros(size(B').*[1,2]);
Btmp(1:2:end) = B.';
B = Btmp.';

% Solve for surface density - this is wrong
tic, X = srskelf_sv_nn(F, B); tsolve = toc;
% X = A\B;
% 
% % Extract part of X corresponding to physical points
% X1 = X(1:2:end);
% X1 = X1./sqrt(area).';
% Y1 = lpcomp_helm_comb_dir(sinfo, zstmp, X1, xyz_out, rank_or_tol);
% 
% % return
% Y_boundary = lpcomp_helm_comb_dir_boundary(sinfo, zstmp, X1, x_or, rank_or_tol);
% % return
% % X2 = X(2:2:end);
% % X2 = X2./sqrt(area).';
% Y2 = lpcomp_helm_comb_dir(sinfo, zdtmp, Y_boundary, xyz_out, rank_or_tol);
% 
% Y = -1j*zk*Y1+Y2;
% % Compare against exact field
% Z = helm_dirichlet_kernel(xyz_out, xyz_in, zstmp, nu2)*q;
% % tmp1 = sqrt(area)'.*X;
% ra =  1; % norm(tmp1);
% e = norm(Z - Y)/ra;
% 
% tmp = helm_dirichlet_kernel(xyz_out, xyz_in, zstmp, nu2);
% % size(tmp)
% 
% fprintf('npts: %d\n', N);
% % fprintf('igeomtype: %d\n',igeomtype);
% fprintf('npatches: %d\n',sinfo.npatches);
% fprintf('norder: %d\n',norder);
% fprintf('zk: %d\n',zk);
% fprintf('time taken for generating quadrature: %d\n',tquad);
% fprintf('time taken for factorization: %d\n',tfac);
% fprintf('time taken for solve: %d\n',tsolve);
fprintf('pde: %10.4e\n',e)

% Examine values
Y(1:2)
Z(1:2)

% Now start scattering test
% diary('off')
% return
% 
% [uinc,xd,xn,thet] = get_uinc(ndir,sinfo,zk);
% tic, Xincsol = srskelf_sv_nn(F,uinc); tsolve = toc;
% 
% exd = exp(-1j*zk*xd);
% dfar = -1j*zpars(3)*zk*xn.*exd;
% sfar = zpars(2)*exd;
% ww = sinfo.wts;
% wwr = repmat(ww,[1,ndir]);
% ufar = (dfar+sfar).*Xincsol.*sqrt(wwr)/4/pi;
% 
% 
% ufar = sum(ufar,1);
% varargout{1} = thet;
% varargout{2} = ufar;
% save(fsol,'ufar','thet','Xincsol');

% edir = norm(Z - Y2)/norm(Z);
% disp(Y2)
% disp(Z)
% fprintf('pde: %10.4e\n',edir)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% end
