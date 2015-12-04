function [lambda, v] = RayleighQuotient(n, A, v, N)
%
% RayleighQuotient(n, A, v, N)
% Approximate an eigenvalue and eigenvector of A using the Rayleigh
% Quotient Iteration.
% n: dimension of A
% A: matrix
% v: approximate eigenvector x, norm(x, inf) = 1
% N: maximum number of iterations
% lambda: approximate eigenvalue

lambda = (v.'*A*v)/(v.'*v);
for k = 1:N
    w = (A - lambda*eye(n))\v;
    v = w/norm(w, inf);
    lambda = (v.'*A*v)/(v.'*v);
end
end
