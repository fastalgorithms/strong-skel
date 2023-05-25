function test_loop()

    % Setup
    addpath('../../fortran_src/')
    addpath('../../../strong-skel/fortran_src/') %git
    addpath('../../')
    addpath('../../sv')
    addpath('../../mv')
    addpath('../../src')
    addpath('../../src/helm_dirichlet/'); %local
    addpath('../../src/helm_sound_hard/'); %local
    run('../../../FLAM/startup.m'); %git

    npu = 10;
    norder = 5;
    occ = 100;
    z_k = 2.0;
    rank_or_tol = 1E-3;
    m = 5;

    % Initialize a wiggly torus
radii = [1.0;2.0;0.25];
scales = [1.0;1.0;1.0];
nnu = npu;
nnv = npu;
nosc = 5;
sinfo = wtorus(radii, scales, nosc, nnu, nnv, norder);

% Pick 'm' points inside object
xyz_in = zeros(3, m);
rng(42);

uu = rand(m,1)*pi/4 -pi/8;
xyz_in(3,:) = radii(1)*sin(uu)*scales(3);
vv = rand(m,1)*2*pi;
rmin = radii(2) + radii(3)*cos(nosc*vv) - radii(1)*abs(cos(uu));
rmax = radii(2) + radii(3)*cos(nosc*vv) + radii(1)*abs(cos(uu));

rr = rand(m,1)*0.1 + 0.45;
rruse = rmin + rr.*(rmax-rmin);
xyz_in(1,:) = rruse.*cos(vv)*scales(1);
xyz_in(2,:) = rruse.*sin(vv)*scales(2);

% Pick 'm' points on parallel surface to the object
xyz_out = zeros(3, m);
uu = rand(m, 1)*2*pi;

zmin = 2.0*radii(1)*scales(3);
zmax = 4.0*radii(1)*scales(3);
xyz_out(3,:) = zmin + rand(m,1)*(zmax-zmin);
rr2 = rand(m,1)*zmax;
xyz_out(1,:) = rr2.*cos(uu);
xyz_out(2,:) = rr2.*sin(uu);

q = rand(m, 1)-0.5 + 1j*(rand(m,1)-0.5);

% Quadrature points
x = sinfo.srcvals(1:3,:);

nu = sinfo.srcvals(10:12,:);
area = sinfo.wts';
N = sinfo.npts;

% Complex parameterizations
zpars = complex([z_k; -1j*z_k; 1.0]);
zstmp = complex([z_k;1.0;0.0]);
zdtmp = complex([z_k;0.0;1.0]);
eps = rank_or_tol;

% Compute the quadrature corrections
tic, Q = helm_neu_near_corr(sinfo, zpars(1:2), eps); tquad=  toc;
S = Q{1};
P = zeros(N,1);
w = whos('S');
fprintf('quad: %10.4e (s) / %6.2f (MB)\n',tquad,w.bytes/1e6)

% Set system matrix, operating on points indices
Afun = @(i, j) Afun_helm_sound_hard(i, j, x, z_k, nu, area, P, S);

% Set proxy function
pxyfun = @(slf, nbr, proxy_dict, l, ctr) pxyfun_helm_sound_hard(x, slf, nbr, proxy_dict, l, ctr, z_k, nu, area);

% Factorize the matrix
opts = struct('verb', 1, 'symm','n', 'z_k', z_k);

F = srskelf_asym_new_pts(Afun, x, occ, rank_or_tol, pxyfun, opts);

end