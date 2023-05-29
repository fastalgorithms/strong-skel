function [nlege] = find_nlege(nlmax,boxsize,dkmax,tol)
    [xlmax,wlmax,ulmax,vlmax] = lege.exps(nlmax);
    sq1 = sqrt(xlmax.^2+(4/5)^2);
    sq2 = sqrt((xlmax+1).^2+(4/5)^2);
    dk  = dkmax*boxsize;
    %f1 = exp(1i*dk*sq1)./(dk*sq1);
    %f2 = exp(1i*dk*sq2)./(dk*sq2);
    f1 = exp(1i*dk*sq1)./(sq1);
    f2 = exp(1i*dk*sq2)./(sq2);
    cfs1 = ulmax*f1;
    cfs2 = ulmax*f2;
    max_coefs = (1+boxsize)*max(abs([cfs1,cfs2].'));
    nlege = find(max_coefs>tol,1,'last');
end

