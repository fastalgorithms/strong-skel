function [varargout] =  ie_fmm3dbie_ss_paper_final_ex(igeomtype,iref,npu,norder,ndir,zk)
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

if(nargin == 0)
    igeomtype = 1;
    npu = 10;
    iref = 0;
    norder = 3;
    zk = 1.0;
    ndir = 100;
end
    
addpath('../fortran_src')
addpath('../')
addpath('../sv')
addpath('../mv')
addpath('../src')
addpath('../src/helm_dirichlet');
run('../../FLAM/startup.m');
occ = 4096;
p = 512;
rank_or_tol = 0.51e-4;
% 
% fname = ['plane-data/diary_igeomtype' int2str(igeomtype) '_iref' int2str(iref) ...
%     '_np' int2str(npu) '_norder' int2str(norder) ...
%     '_ndir' int2str(ndir) '_zk' num2str(zk) '_ss.dat'];
% fsol = ['plane-data/sol_igeomtype' int2str(igeomtype) '_iref' int2str(iref) ...
%     '_np' int2str(npu) '_norder' int2str(norder) ...
%     '_ndir' int2str(ndir) '_zk' num2str(zk) '_ss.dat'];
% 
% diary(fname);

if(igeomtype == 1)


    % Initialize the discretization points on the sphere and the proxy 
    % sphere.  We use random points on the proxy sphere for simplicity.
    radii = [1.0;2.0;0.25];
    scales = [1.2;1.0;1.7];

    nnu = npu;
    nnv = npu;
    nosc = 5;
    sinfo = wtorus(radii,scales,nosc,nnu,nnv,norder);
    
    m = 50;
    rng(42);
    uu = rand(m,1)*2*pi;
    vv = rand(m,1)*2*pi;
    rr = rand(m,1)*0.67;
    xyz_in = zeros(3,m);
    xyz_in(1,:) = (rr.*cos(uu) + 2 + 0.25*cos(5*vv)).*cos(vv)*1.2;
    xyz_in(2,:) = (rr.*cos(uu) + 2 + 0.25*cos(5*vv)).*sin(vv)*1.0;
    xyz_in(3,:) = rr.*sin(uu)*1.7;
    
    xyz_out = zeros(3,m);
    uu = rand(m,1)*2*pi;
    vv = rand(m,1)*2*pi;
    rr = rand(m,1)*0.67 + 1.33;

    xyz_out(1,:) = (rr.*cos(uu) + 2 + 0.25*cos(5*vv)).*cos(vv)*1.2;
    xyz_out(2,:) = (rr.*cos(uu) + 2 + 0.25*cos(5*vv)).*sin(vv)*1.0;
    xyz_out(3,:) = rr.*sin(uu)*1.7;


else
   sinfo = plane_a380(iref); 
   xyz_in = load('fuselage-sources_matlab.txt');
   xyz_in = xyz_in.';
   
   tmp = load('targ_plane_ext.dat');
   tmp = tmp.';
   xyz_out = tmp(1:3,:);
end
x = sinfo.srcvals(1:3,:);
nu = sinfo.srcvals(10:12,:);
area = sinfo.wts';
N = sinfo.npts;

proxy = randn(3,p);
proxy = 1.5*bsxfun(@rdivide,proxy,sqrt(sum(proxy.^2)));

% Compute the quadrature corrections

zpars = complex([zk; -1j*zk; 1.0]);
zstmp = complex([zk;1.0;0.0]);
zdtmp = complex([zk;0.0;1.0]);
eps = rank_or_tol;
tic, S = helm_near_corr(sinfo,zpars,eps); tquad=  toc;
P = zeros(N,1);
w = whos('S');
fprintf('quad: %10.4e (s) / %6.2f (MB)\n',tquad,w.bytes/1e6)


% Factor the matrix using skeletonization (verbose mode)
opts = struct('verb',1,'symm','n','zk',zk);

Afun_use = @(i,j) Afun_helm_dirichlet(i,j,x,zpars,nu,area,P,S);
pxyfun_use = @(x,slf,nbr,proxy,l,ctr) pxyfun_helm_dirichlet(x,slf,nbr,proxy,l,ctr,zpars,nu,area);
tic, F = srskelf_asym_new(Afun_use,x,occ,rank_or_tol,pxyfun_use,opts); tfac = toc;
w = whos('F');
fprintf([repmat('-',1,80) '\n'])
fprintf('mem: %6.4f (GB)\n',w.bytes/1048576/1024)
%save(fname2,'F');

% 
q = rand(m,1)-0.5; + 1j*(rand(m,1)-0.5);
nu2 = zeros(3,m);
B = helm_dirichlet_kernel(x,xyz_in,zstmp,nu2)*q;
B = B.*sqrt(area).';

% Solve for surface density
tic, X = srskelf_sv_nn(F,B); tsolve = toc;
X = X./sqrt(area).';
Y = lpcomp_helm_comb_dir(sinfo,zpars,X,xyz_out,rank_or_tol);

% Compare against exact field
Z = helm_dirichlet_kernel(xyz_out,xyz_in,zstmp,nu2)*q;
tmp1 = sqrt(area)'.*X;
ra = norm(tmp1);
e = norm(Z - Y)/ra;

fprintf('npts: %d\n',N);
fprintf('igeomtype: %d\n',igeomtype);
fprintf('npatches: %d\n',sinfo.npatches);
fprintf('norder: %d\n',norder);
fprintf('zk: %d\n',zk);
fprintf('time taken for generating quadratue: %d\n',tquad);
fprintf('time taken for factorization: %d\n',tfac);
fprintf('time taken for solve: %d\n',tsolve);
fprintf('pde: %10.4e\n',e)


% Now start scattering test


diary('off')
return

[uinc,xd,xn,thet] = get_uinc(ndir,sinfo,zk);
tic, Xincsol = srskelf_sv_nn(F,uinc); tsolve = toc;

exd = exp(-1j*zk*xd);
dfar = -1j*zpars(3)*zk*xn.*exd;
sfar = zpars(2)*exd;
ww = sinfo.wts;
wwr = repmat(ww,[1,ndir]);
ufar = (dfar+sfar).*Xincsol.*sqrt(wwr)/4/pi;


ufar = sum(ufar,1);
varargout{1} = thet;
varargout{2} = ufar;
save(fsol,'ufar','thet','Xincsol');

% edir = norm(Z - Y2)/norm(Z);
% disp(Y2)
% disp(Z)
% fprintf('pde: %10.4e\n',edir)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end

