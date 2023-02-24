function [varargout] =  test2(npu, norder, occ, zk, rank_or_tol, m)

%%%%% a basic test code for constructing operators for transmission 
%%%%% problems. Not supposed to be used more generally. Compares with 
%%%%% analytic solution (exterior only...)

addpath('../../fortran_src/')
addpath('../../../strong-skel/fortran_src/') %git
addpath('../../')
addpath('../../sv')
addpath('../../mv')
addpath('../../src')
addpath('../../src/helm_dirichlet/'); %local
addpath('../../src/helm_sound_hard/'); %local
addpath('../../src/helm_transmission/'); %local
run('../../../FLAM/startup.m'); %git

eps1 = 2.1;
eps2 = 1.4;
epss = [eps1,eps2];

mu1 = 2;
mu2 = 3;
mus = [mu1,mu2];

if(nargin == 0)

    npu = 10;
    norder = 5;
    occ = 1000;
    zk = 1.0;
    rank_or_tol = 5e-7;
    rank_or_tol = 1E-6;
    m = 10;
    m2= 12;
end

a = 1.2;
b = 1;
c = 1.4;
ifc = 0;
rmax = min(0.5,pi/(abs(zk)*max(sqrt(epss.*mus))));
sinfo = ellipsoid(a,b,c,rmax,ifc,norder);

% Pick 'm' points inside object
xyz_in = zeros(3, m);
ifread = 0;
ifwrite = 0;
if(~ifread)
    rng(42);

    xyz_in = randn(3,m);
    xyz_in = (1-rmax)*xyz_in.*(rand(1,m)./vecnorm(xyz_in,2));
    xyz_in = diag([a,b,c])*xyz_in;
% Pick 'm' points outside surface to the object
    xyz_out = randn(3, m2);
    rads = (0.5*rand(1,m2)+1.5);
    xyz_out = xyz_out.*(rads./vecnorm(xyz_out,2));
    xyz_out = diag([a,b,c])*xyz_out;
    
    q = rand(m, 1)-0.5 + 1j*(rand(m,1)-0.5);
    %q(1) = 0;
    %q=rand(1)+1j*rand(1);
else
    xyz_in = readmatrix('in.txt').';
    xyz_out = readmatrix('out.txt').';
    qq = readmatrix('q.txt');
    q = qq(:,1) + 1j*qq(:,2);

end

if(ifwrite)
   writematrix(xyz_in.','in.txt');
   writematrix(xyz_out.','out.txt');
   writematrix([real(q) imag(q)],'q.txt');
end


% Quadrature points
x = sinfo.srcvals(1:3,:);
x = repelem(x, 1,2);
x_or = sinfo.srcvals(1:3,:);

% figure(1)
% clf
% scatter3(x_or(1,:),x_or(2,:),x_or(3,:),'k.'); hold on;
% scatter3(xyz_in(1,:),xyz_in(2,:),xyz_in(3,:),'ro');

nu = sinfo.srcvals(10:12,:);
area = sinfo.wts';
N = sinfo.npts;

% Complex parameterizations
zpars = complex([zk;eps1;mu1;eps2;mu2]);

% not sure if necessary...
zstmp = complex([zk;1.0;0.0]);
zdtmp = complex([zk;0.0;1.0]);
eps = rank_or_tol;

% Compute the quadrature corrections
tic, S = helm_trans_near_corr(sinfo,zpars,eps); tquad =toc;
P = zeros(N,1);
w = whos('S');
fprintf('quad: %10.4e (s) / %6.2f (MB)\n',tquad,w.bytes/1e6)

% Set system matrix
Afun_use = @(i, j) Afun_helm_transmission_wrapper(i, j, x_or, zpars, nu, area, P, S);

% Set proxy function
% pxyfun_use = @(xx, slf, nbr, proxy, l, ctr) pxyfun_helm_sound_hard_wrapper(xx, slf, nbr, proxy, l, ctr, x_or, zk, nu, area);
pxyfun_use = [];
% Factorize the matrix
opts = struct('verb', 1, 'symm','n', 'zk', zk);

%tic, F = srskelf_asym_new(Afun_use, x, occ, rank_or_tol, pxyfun_use, opts); tfac = toc;
%w = whos('F');
%fprintf([repmat('-',1,80) '\n'])
%fprintf('mem: %6.4f (GB)\n',w.bytes/1048576/1024)


zk0 = zpars(1)*sqrt(zpars(4)*zpars(5));
deps0 = zpars(4);
zk  = zpars(1)*sqrt(zpars(2)*zpars(3));
deps  = zpars(2);

[u0,du0] = helm_transmission_kernel(x_or, xyz_in, zk, nu);
B1 =  (u0*q).*sqrt(area).';
B2 = -(du0*q/deps).*sqrt(area).';
Btmp = zeros(size(B1').*[1,2]);
Btmp(1:2:end) = B1;
Btmp(2:2:end) = B2;
B = Btmp.';

%B(1:2:end) = sqrt(area);
%B(2:2:end) = 0;

% Solve for surface density
%tic, X = srskelf_sv_nn(F, B); tsolve = toc;

 A = Afun_use(1:2*N, 1:2*N);
 X_test = A \ B;

 dens1 = X_test(1:2:end)./sqrt(area.');
 dens2 = X_test(2:2:end)./sqrt(area.');
 
 [G,dG] = helm_transmission_kernel(x_or, xyz_out, zk, nu);
 u0comp  =  sum(G.*dens2.*(area.'),1);
 u0dcomp =  sum(dG.*dens1.*(area.'),1);
 u0 = deps^2*u0comp + deps*u0dcomp;
 
 nvecs = zeros(size(xyz_in));
[u0true,du0true] = helm_transmission_kernel(xyz_in, xyz_out, zk, nvecs);
 u0true = (u0true.')*q;
 
end