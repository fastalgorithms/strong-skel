
rng('default');
rng(1);

% Number of points
N = 100;
% Block dimensions
m = 2;
n = 2;
x = randn(N, 2);

A = LaplaceBlock(x, m, n);

tol = 1e-4;

% [sk, rd, T, niter] = id(A, tol);
%
% a = A(:, rd);
% b = A(:,sk)*T + tol;

% e = norm(a-b, "fro");

% Reshape A to have blocks of dim ((m*n)x1) instead of (mxn)
Areshaped = zeros(m*n*N, N);

for i = 1:N
    for j = 1:N
        tmp_block = A(m*(i-1)+1:m*i, n*(j-1)+1:n*j);
        reshaped_block = reshape(tmp_block, 1, [])';
        Areshaped(m*n*(i-1)+1:m*n*i, (j-1)+1:j) = reshaped_block;
    end
end

% Try computing ID on reshaped A

[sk, rd, T] = id(Areshaped, tol);

n_rd_or = length(rd)*n;
n_sk_or = length(sk)*n;
rd_or = zeros(1, n_rd_or);
sk_or = zeros(1, n_sk_or);

% Find sk, rd and T in terms of redundant indices
for i = 1:length(sk)
    tmp = sk(i);
    for j = 1:n
        sk_or(n*(i-1)+j) = (n*tmp)+j-n;
    end
end


for i = 1:length(rd)
    tmp = rd(i);
    for j = 1:n
        rd_or(n*(i-1)+j) = (n*tmp)+j-n;
    end
end

T_or = zeros(size(T, 1)*n, size(T, 2)*m);
for i = 1:size(T, 1)
    for j = 1:size(T, 2)
        T_or((i-1)*n+1:i*n, (j-1)*m+1:j*m) = repmat(T(i, j), n, m);
    end
end

% T
% T_or_2 = repelem(T, 2, 2)
% T_or
% size(sk)
% size(sk_or)
% size(T)
% size(T_or)
% size(A(:, sk_or))
% Try and reconstruct block matrix from ID over point indices

a = A(:, rd_or);
b = A(:, sk_or)*T_or + tol;
a(1:10, 1)'
1/n * b(1:10, 1)'

% repmat(T(1, 1), M, M)

function [A] = LaplaceBlock(x, m, n)
% Some rank structured block matrix for testing
    N = length(x);
    A = zeros(m*N, n*N);

    for i=1:length(x)
        for j=1:length(x)

            % Laplace kernel
            tmp = 1./sqrt(sum(abs((x(i)-x(j)))^2));

            % Block value
            tmp_block = zeros(m, n);
            for k = 1:m
                for l = 1:n
                    tmp_block(k, l) = tmp;
                    if isinf(tmp)
                        tmp_block(k, l) = 0.;
                    end
                end
            end

            A(m*(i-1)+1:m*i, n*(j-1)+1:n*j) = tmp_block;
        end
    end
end

