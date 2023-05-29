function [proxy_dict] = init_proxy_dict(boxsize,opts,tol)

    nterms = h3dterms(boxsize,opts.zk,tol);
    dkmax = max(abs(opts.zk));
    RR = 5/2;
    nleg = 2*ceil((nterms+1));
    nleg_fmm = nleg
    nlmax  = 40 + ceil(RR*boxsize*dkmax/pi);
    nlmax
    
    ifdone = false;
    while (~ifdone)
        [nleg] = find_nlege(nlmax,boxsize,dkmax,tol)
        ifdone = nleg<nlmax;
        nlmax = nlmax*2;
    end
    nleg = nlmax;
    [xleg,wleg] = lege.exps(nleg);

    
    [XL,YL] = meshgrid(xleg,xleg);
    [WX,WY] = meshgrid(wleg,wleg);
    XL = RR*XL(:);
    YL = RR*YL(:);
    ZL = ones(size(XL));
    WL = RR*sqrt(WX(:).*WY(:));
    
    DO = ones(size(ZL));
    DZ = zeros(size(ZL));
    
    %
    ppts = [XL,YL,RR*ZL];
    nrml = [DZ,DZ,DO];
    wpts = WL;
    %
    ppts = [ppts; XL,YL,-RR*ZL];
    nrml = [nrml;DZ,DZ,-DO];
    wpts = [wpts;WL];
    %
    ppts = [ppts; XL,-RR*ZL,YL];
    nrml = [nrml;DZ,-DO,DZ];
    wpts = [wpts;WL];
    %
    ppts = [ppts; XL,RR*ZL,YL];
    nrml = [nrml;DZ,DO,DZ];
    wpts = [wpts;WL];
    %
    ppts = [ppts; -RR*ZL,XL,YL];
    nrml = [nrml;-DO,DZ,DZ];
    wpts = [wpts;WL];
    %
    ppts = [ppts; RR*ZL,XL,YL];
    nrml = [nrml;DO,DZ,DZ];
    wpts = [wpts;WL];
    
    proxy_dict = [];
    proxy_dict.proxy = ppts';
    proxy_dict.weigt = wpts';
    proxy_dict.norms = nrml';
end

