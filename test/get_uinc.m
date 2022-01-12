function [uinc,xd,xn,thet] = get_uinc(ndir,sinfo,zk)
    x = sinfo.srcvals(1:3,:);
    thet = 2*pi/ndir:2*pi/ndir:2*pi;

    thet = thet(:);
    ct = cos(thet);
    st = sin(thet);
    xx = x(1,:);
    xx = xx(:);
    yy = x(2,:);
    yy = yy(:);
   
    rnx = sinfo.srcvals(10,:);
    rnx = rnx(:);
    rny = sinfo.srcvals(11,:);
    rny = rny(:);
    xd = xx*ct' + yy*st';
    xn = rnx*ct' + rny*ct';

    uinc = - exp(1j*zk*xd);
    ww = sinfo.wts;
    wwr = repmat(ww,[1,ndir]);
    uinc = uinc.*sqrt(wwr);

end