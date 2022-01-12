% Test plane wave computation and far-field computation

sinfo  = plane_a380(0);
x = sinfo.srcvals(1:3,:);
ndir = 100;
thet = 2*pi/ndir:2*pi/ndir:2*pi;

thet = thet(:);
ct = cos(thet);
st = sin(thet);
xx = x(1,:);
xx = xx(:);
yy = x(2,:);
yy = yy(:);
zz = x(3,:);
zz = zz(:);

xd = xx*ct' + yy*st';

zk = 1.0;
uinc = - exp(1j*zk*xd);
ww = sinfo.wts;
wwr = repmat(ww,[1,ndir]);
uinc = uinc.*sqrt(wwr);

