function test1(ik,npu,norder)
radii = [1.0;2.0;0.25];
scales = [1.2;1.0;1.7];
nu = npu;
nv = npu;
nosc = 5;
norder = norder;
if(ik == 1)
    zk = 0.97;
else
    zk = 0.97*nu/10;
end
fname = ['diary_ik' int2str(ik) '_np' int2str(npu) '_norder' int2str(norder) '_qtime.dat'];
diary(fname);
S = wtorus(radii,scales,nosc,nu,nv,norder);
zpars = complex([zk; 1.0; 0.0]);
eps = 0.51e-7;
tic, spmat = helm_near_corr(S,zpars,eps); toc;
zpars = complex([zk; 0.0; 1.0]);
tic, spmat = helm_near_corr(S,zpars,eps); toc;
zpars = complex([zk; 1j*zk; 1.0]);
tic, spmat = helm_near_corr(S,zpars,eps); toc;
zpars = complex([zk; -1j*zk; 1.0]);
tic, spmat = helm_near_corr(S,zpars,eps); toc;
diary('off');
exit;
end
