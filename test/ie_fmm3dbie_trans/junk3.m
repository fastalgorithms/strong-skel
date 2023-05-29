ns =1:400;
xs = 0:0.01:(2*pi);

[N,X] = meshgrid(ns,xs);

cfs = (-1).^N./N.^(0.6).*exp(1i*X.*N);

vals = sum(cfs,2);
plot(xs,imag(vals))