function ie_fmm3dbie_sqrtscaling_stor_load()
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
%                 used to automatically select the number of skeletons.
addpath('../fortran_src')
addpath('../')
addpath('../sv')
addpath('../mv')
run('../../FLAM/startup.m');
occ = 256;
p = 512;
rank_or_tol = 0.51e-6;
maxNumCompThreads(1);

zk = 0.97;
npu = 20;
norder = 7;

ifload = 1;
ifstor = 0;

fnameq = ['quad_wtorus_np' int2str(npu) '_norder' int2str(norder) '.mat'];
fnamef = ['factorization_wtorus_np' int2str(npu) '_norder' int2str(norder) '.mat'];

% Initialize the discretization points on the sphere and the proxy 
% sphere.  We use random points on the proxy sphere for simplicity.
radii = [1.0;2.0;0.25];
scales = [1.2;1.0;1.7];

nnu = npu;
nnv = npu;
nosc = 5;
sinfo = wtorus(radii,scales,nosc,nnu,nnv,norder);
x = sinfo.srcvals(1:3,:);
nu = sinfo.srcvals(10:12,:);
area = sinfo.wts';
N = sinfo.npts;

proxy = randn(3,p);
proxy = 1.5*bsxfun(@rdivide,proxy,sqrt(sum(proxy.^2)));

zpars = complex([zk; -1j*zk; 1.0]);
zstmp = complex([zk;1.0;0.0]);
zdtmp = complex([zk;0.0;1.0]);
eps = 0.51e-6;

fprintf('Starting quadrature generation:\n');
if (ifload == 0)
   tic, S = helm_near_corr(sinfo,zpars,eps); tquad=  toc;
else
   tic, load(fnameq,'S'); tquad=toc;
end

if (ifstor == 1)
   save(fnameq,'S','-v7.3');
end
P = zeros(N,1);
w = whos('S');
fprintf('quad: %10.4e (s) / %6.2f (MB)\n',tquad,w.bytes/1e6)


% Factor the matrix using skeletonization (verbose mode)
opts = struct('verb',1,'symm','n','zk',zk);

if(ifload == 0)
   tic, F = srskelf_asym_new(@Afun,x,occ,rank_or_tol,@pxyfun,opts); tfac = toc;
else
   tic, load(fnamef,'F'); tfac = toc;
end

if(ifstor == 1)
   save(fnamef,'F','-v7.3');
end
w = whos('F');
fprintf([repmat('-',1,80) '\n'])
fprintf('mem: %6.4f (GB)\n',w.bytes/1048576/1024)

m = 50;
rng(42);
uu = rand(m,1)*2*pi;
vv = rand(m,1)*2*pi;
rr = rand(m,1)*0.67;
xyz_in = zeros(3,m);
xyz_in(1,:) = (rr.*cos(uu) + 2 + 0.25*cos(5*vv)).*cos(vv)*1.2;
xyz_in(2,:) = (rr.*cos(uu) + 2 + 0.25*cos(5*vv)).*sin(vv)*1.0;
xyz_in(3,:) = rr.*sin(uu)*1.7;
% 
% m = 2;
% src = [0.11,0.13;-2.13,2.1;0.05,-0.01];
q = rand(m,1)-0.5; + 1j*(rand(m,1)-0.5);
%q = [1.0;1.0+1.0*1j];
nu2 = zeros(3,m);
B = Kfun(x,xyz_in,zstmp,nu2)*q;
B = B.*sqrt(area).';

% Solve for surface density
tic, X = srskelf_sv_nn(F,B); tsolve = toc;
X = X./sqrt(area).';

% A2 = Afun(1:N,1:N);
% X2 = A2\B;

% Evaluate field at interior targets
xyz_out = zeros(3,m);
uu = rand(m,1)*2*pi;
vv = rand(m,1)*2*pi;
rr = rand(m,1)*0.67 + 1.33;

xyz_out(1,:) = (rr.*cos(uu) + 2 + 0.25*cos(5*vv)).*cos(vv)*1.2;
xyz_out(2,:) = (rr.*cos(uu) + 2 + 0.25*cos(5*vv)).*sin(vv)*1.0;
xyz_out(3,:) = rr.*sin(uu)*1.7;


%trg = [31.17,6.13;-0.03,-4.1;3.15,22.2];
%Y = bsxfun(@times,Kfun(trg,x,zpars,nu),area)*X;
Y = lpcomp_helm_comb_dir(sinfo,zpars,X,xyz_out,rank_or_tol);
%Y2 = bsxfun(@times,Kfun(trg,x,zpars,nu),area)*X2;

% Compare against exact field
Z = Kfun(xyz_out,xyz_in,zstmp,nu2)*q;
tmp1 = sqrt(area)'.*X;
ra = norm(tmp1);
e = norm(Z - Y)/ra;

fprintf('npts: %d\n',N);
fprintf('npatches: %d\n',sinfo.npatches);
fprintf('norder: %d\n',norder);
fprintf('zk: %d\n',zk);
fprintf('time taken for generating quadratue: %d\n',tquad);
fprintf('time taken for factorization: %d\n',tfac);
fprintf('time taken for solve: %d\n',tsolve);
fprintf('pde: %10.4e\n',e)

% edir = norm(Z - Y2)/norm(Z);
% disp(Y2)
% disp(Z)
% fprintf('pde: %10.4e\n',edir)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function K = Kfun(x,y,zpuse,nuuse)
% KFUN(X,Y,zpars,NU) computes the Helmholtz potential evaluated
% pairwise between points in X and points in Y (does not handle the
% singularity).  zpars says which layer potential to use
% using the surface normal vectors in NU.
dx = bsxfun(@minus,x(1,:)',y(1,:));
dy = bsxfun(@minus,x(2,:)',y(2,:));
dz = bsxfun(@minus,x(3,:)',y(3,:));
dr = sqrt(dx.^2 + dy.^2 + dz.^2);
zexp = exp(1j*zpuse(1)*dr);
slp = 1/(4*pi)./dr.*zexp;
rdotn = bsxfun(@times,dx,nuuse(1,:)) + bsxfun(@times,dy,nuuse(2,:)) + ...
          bsxfun(@times,dz,nuuse(3,:));
dlp = 1/(4*pi).*rdotn./dr.^3.*zexp.*(1.0-1j*zpuse(1)*dr);
K = zpuse(2)*slp + zpuse(3)*dlp;
K(dr == 0) = 0;
end

function A = Afun(i,j)
% AFUN(I,J) computes entries of the matrix A to be factorized at the
% index sets I and J.  This handles the near-field correction.
if isempty(i) || isempty(j)
  A = zeros(length(i),length(j));
  return
end
[I,J] = ndgrid(i,j);
A = bsxfun(@times,Kfun(x(:,i),x(:,j),zpars,nu(:,j)),area(j));
M = spget(i,j);
idx = abs(M) ~= 0;
A(idx) = M(idx);
A(I == J) = A(I == J) + 0.5*zpars(3);

A = bsxfun(@times,sqrt(area(i)).',A);
A = bsxfun(@times,A,1.0./sqrt(area(j)));
end

% proxy function
function [Kpxy,nbr] = pxyfun(x,slf,nbr,proxy,l,ctr)
% PXYFUN(X,SLF,NBR,L,CTR) computes interactions between the points
% X(:,SLF) and the set of proxy points by scaling the proxy sphere to 
% appropriately contain a box at level L centered at CTR and then
% calling KFUN
pxy = bsxfun(@plus,proxy*l,ctr');
Kpxy1 = Kfun(pxy,x(:,slf),zstmp,nu(:,slf));
Kpxy1 = bsxfun(@times,Kpxy1,sqrt(area(slf)));
Kpxy3 = bsxfun(@times,Kfun(pxy,x(:,slf),zpars,nu(:,slf)),sqrt(area(slf)));
Kpxy = [Kpxy1;Kpxy3];
dx = x(1,nbr) - ctr(1);
dy = x(2,nbr) - ctr(2);
dz = x(3,nbr) - ctr(3);
dist = sqrt(dx.^2 + dy.^2 + dz.^2);
nbr = nbr(dist/l < 1.5);
end

function A = spget(I_,J_)
% SPGET(I_,J_) computes entries of a sparse matrix of near-field
% corrections that should be added to the kernel matrix, as used in
% AFUN.
m_ = length(I_);
n_ = length(J_);
[I_sort,E] = sort(I_);
P(I_sort) = E;
A = zeros(m_,n_);
[I_,J_,S_] = find(S(:,J_));
idx = ismembc(I_,I_sort);
I_ = I_(idx);
J_ = J_(idx);
S_ = S_(idx);
A(P(I_) + (J_ - 1)*m_) = S_;
end
end

