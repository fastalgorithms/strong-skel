
function A = Afun_helm_dirichlet(i,j,x,zpars,nu,area,P,S)
% AFUN(I,J) computes entries of the matrix A to be factorized at the
% index sets I and J.  This handles the near-field correction.
if isempty(i) || isempty(j)
  A = zeros(length(i),length(j));
  return
end
[I,J] = ndgrid(i,j);
A = bsxfun(@times,helm_dirichlet_kernel(x(:,i),x(:,j),zpars,nu(:,j)),area(j));
M = spget_quadcorr(i,j,P,S);
idx = abs(M) ~= 0;
A(idx) = M(idx);
A(I == J) = A(I == J) + 0.5*zpars(3);

A = bsxfun(@times,sqrt(area(i)).',A);
A = bsxfun(@times,A,1.0./sqrt(area(j)));
end
