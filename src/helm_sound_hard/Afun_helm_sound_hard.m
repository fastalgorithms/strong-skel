
function A = Afun_helm_sound_hard(i, j, x, zpars, nu, area, P, S)
% AFUN(I,J) computes entries of the matrix A to be factorized at the
% index sets I and J.  This handles the near-field correction.

if isempty(i) || isempty(j)
  A = zeros(2*length(i), 2*length(j));
  return
end

% Complex wavenumber
z_k = zpars(1);

% Quadrature corrections
% M = spget_quadcorr(i, j, P, S);
% idx = abs(M) ~= 0;
% A(idx) = M(idx);

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
         A22 = complex(-1);
         
         A11(ii == jj) = A11(ii == jj) + (1j*z_k*0.5 - 0.25);

         block = [A11 A12; A21 A22];
       
         r_idx = ii*2 - 1;
         c_idx = jj*2 - 1;
         A(r_idx:r_idx+1, c_idx:c_idx+1) = block;

    end
end

% A = bsxfun(@times,sqrt(area(i)).',A);
% A = bsxfun(@times,A,1.0./sqrt(area(j)));
end
