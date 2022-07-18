
function A = Afun_helm_sound_hard(i, j, x, zpars, nu, area, P, S)
% AFUN(I,J) computes entries of the matrix A to be factorized at the
% index sets I and J.  This handles the near-field correction.

if isempty(i) || isempty(j)
  A = zeros(2*length(i), 2*length(j));
  return
end

% Complex wavenumber
z_k = zpars(1);

% Assemble (block) system matrix for sound hard problem
A = zeros(2*length(i), 2*length(j));
for ii = 1:length(i)
   
    for jj = 1:length(j)
        
         targets = x(:, ii);
         nx = nu(:, ii);
         sources = x(:, jj);
         K_prime = bsxfun(@times, helm_sound_hard_kernel(targets, sources, z_k, nx), area(jj));

         A11 =  - 1j*z_k*K_prime;
         A12 = K_prime;
         A21 = K_prime;
         A22 = 0;
         
         % Compute quadrature corrections
         M = spget_quadcorr(ii, jj, P, S);
         idx = abs(M) ~= 0;

         % Add quadrature corrections to each block
         A11(idx) = M(idx);
         A12(idx) = M(idx);
         A21(idx) = M(idx);
         A22(idx) = M(idx);

         A11(ii == jj) = A11(ii == jj) + (1j*z_k*0.5 - 0.25);
         A22(ii == jj) = complex(-1);

         % Compute weights for each block
         A11 = bsxfun(@times, sqrt(area(ii)).',A11);
         A11 = bsxfun(@times, A11,1.0./sqrt(area(jj)));
         A12 = bsxfun(@times, sqrt(area(ii)).',A12);
         A12 = bsxfun(@times, A12, 1.0./sqrt(area(jj)));
         A21 = bsxfun(@times, sqrt(area(ii)).',A21);
         A21 = bsxfun(@times, A21, 1.0./sqrt(area(jj)));
         A22 = bsxfun(@times, sqrt(area(ii)).',A22);
         A22 = bsxfun(@times, A22, 1.0./sqrt(area(jj)));

         block = [A11 A12; A21 A22];

         r_idx = ii*2 - 1;
         c_idx = jj*2 - 1;
         A(r_idx:r_idx+1, c_idx:c_idx+1) = block;

    end
end


end
