function [varargout] =  ie_fmm3dbie_sound_hard(npu, norder, occ, zk, rank_or_tol, m)
% IE_SOUND_HARD  An example usage of strong skeletonization, solving a
%  second-kind integral equation (Helmholtz combined field potential) on a
%  wiggly torus for the sound hard scattering problem.
%
%  - NPU:         The number of patches in the u direction of the torus, defines
%                 the discretization level of the torus
%
%  - NORDER:      The order of expansion on each patch.
%
%
%  - OCC:         The occupancy parameter, specifying the maximum number of
%                 points a node in the octree can contain before it is
%                 subdivided.  This therefore gives an upper bound on the
%                 number of points in a leaf node.
%
%  - ZK:          The complex wavenumber.
%
%  - RANK_OR_TOL: If a natural number, the maximum number of skeletons to
%                 select during a single box of skeletonization.  If a
%                 float between 0 and 1, an approximate relative tolerance
%                 used to automatically select the number of skeletons.
%
%  - m:           The number of test points at which to compare the computed
%                 solution with the exact solution, calculated via potential
%                 theory.

addpath('../../fortran_src')
addpath('../../')
addpath('../../sv')
addpath('../../mv')
addpath('../../src')
addpath('../../src/helm_dirichlet/');
addpath('../../src/helm_sound_hard/');
run('../../../FLAM/startup.m');

if(nargin == 0)
    npu = 10;
    norder = 3;
    occ = 50;
    zk = 1.0;
    rank_or_tol = 5e-7;
    m = 10;
end

% Seed for random numbers
rng(42);

% Initialize a wiggly torus
radii = [1.0;2.0;0.25];
scales = [1.2;1.0;0.5];
nnu = npu;
nnv = npu;
nosc = 5;
sinfo = wtorus(radii, scales, nosc, nnu, nnv, norder);

% Pick 'm' points on object surface
xyz_in = zeros(3, m);
uu = rand(m,1)*2*pi;
vv = rand(m,1)*2*pi;
rr = rand(m,1)*0.67;

xyz_in(1,:) = (rr.*cos(uu) + 2 + 0.25*cos(5*vv)).*cos(vv)*1.2;
xyz_in(2,:) = (rr.*cos(uu) + 2 + 0.25*cos(5*vv)).*sin(vv)*1.0;
xyz_in(3,:) = rr.*sin(uu)*1.7;

% Pick 'm' points on parallel surface to the object
xyz_out = zeros(3, m);
uu = rand(m, 1)*2*pi;
vv = rand(m, 1)*2*pi;
rr = rand(m, 1)*0.67 + 10.33;

xyz_out(1,:) = (rr.*cos(uu) + 2 + 0.25*cos(5*vv)).*cos(vv)*1.2;
xyz_out(2,:) = (rr.*cos(uu) + 2 + 0.25*cos(5*vv)).*sin(vv)*1.0;
xyz_out(3,:) = rr.*sin(uu)*1.7;

% Quadrature points
x = sinfo.srcvals(1:3,:);
x = repmat(x, [1,2]);
x_or = sinfo.srcvals(1:3,:);

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
Afun_use = @(i, j) Afun_helm_sound_hard_wrapper(i, j, x_or, zk, nu, area, P, S);

% Set proxy function
pxyfun_use = @(x, slf, nbr, proxy, l, ctr) pxyfun_helm_sound_hard_wrapper(x_or, slf, nbr, proxy, l, ctr, zk, nu, area);

% Factorize the matrix
opts = struct('verb', 1, 'symm','n', 'zk', zk);
tic, F = srskelf_asym_new(Afun_use, x, occ, rank_or_tol, pxyfun_use, opts); tfac = toc;
w = whos('F');
fprintf([repmat('-',1,80) '\n'])
fprintf('mem: %6.4f (GB)\n',w.bytes/1048576/1024)

% Compute Neumann data - g = \Del_x S(x,y)
q = rand(m, 1)-0.5 + 1j*(rand(m,1)-0.5);
B = helm_sound_hard_kernel(x_or, xyz_in, zk, nu)*q;
B = B.*sqrt(area).';

% Add zero data due to system - [g, 0, ..., g, 0]
Btmp = zeros(size(B').*[1,2]);
Btmp(1:2:end) = B.';
B = Btmp.';

% Solve for surface density
tic, X = srskelf_sv_nn(F, B); tsolve = toc;

% Extract part of X corresponding to physical points
X1 = X(1:2:end);
X1 = X1./sqrt(area).';

% Compute potential using combined representation
Y1 = lpcomp_helm_comb_dir(sinfo, zstmp, X1, xyz_out, rank_or_tol);
Y_boundary = lpcomp_helm_comb_dir_boundary(sinfo, zstmp, X1, x_or, rank_or_tol);
Y2 = lpcomp_helm_comb_dir(sinfo, zdtmp, Y_boundary, xyz_out, rank_or_tol);
Y = -1j*zk*Y1+Y2;

% Compare against exact potential
nu2 = zeros(3, m);
Z = helm_dirichlet_kernel(xyz_out, xyz_in, zstmp, nu2)*q;

% Compute a relative error at 'm' points
tmp1 = sqrt(area)'.*X1;
ra = norm(tmp1);
e = norm(Z - Y);

fprintf('npts: %d\n', N);
fprintf('npatches: %d\n',sinfo.npatches);
fprintf('norder: %d\n',norder);
fprintf('zk: %d\n',zk);
fprintf('time taken for generating quadrature: %d\n',tquad);
fprintf('time taken for factorization: %d\n',tfac);
fprintf('time taken for solve: %d\n',tsolve);
fprintf('pde: %10.4e\n',e);

% Now start scattering test
% [uinc, xd, xn, thet] = get_uinc(ndir, sinfo, zk);
% tic, Xincsol = srskelf_sv_nn(F,uinc); tsolve = toc;

% exd = exp(-1j*zk*xd);
% dfar = -1j*zpars(3)*zk*xn.*exd;
% sfar = zpars(2)*exd;
% ww = sinfo.wts;
% wwr = repmat(ww,[1,ndir]);
% ufar = (dfar+sfar).*Xincsol.*sqrt(wwr)/4/pi;

% ufar = sum(ufar,1);
% varargout{1} = thet;
% varargout{2} = ufar;
% save(fsol,'ufar','thet','Xincsol');

% edir = norm(Z - Y2)/norm(Z);
% disp(Y2)
% disp(Z)
% fprintf('pde: %10.4e\n',edir)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%quit;
quit;
end
