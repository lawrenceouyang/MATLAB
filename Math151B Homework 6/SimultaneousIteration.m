function [lambda, W] = SimultaneousIteration(A, V, N)
%
% SimultaneousIteration(A, V, N)
% Approximate an eigenvalue and eigenvector of A using the Simultaneous
% Iteration method.
% A: matrix
% V: matrix of orthogonal column entries
% N: maximum number of iterations
% lambda: approximate eigenvalue vector
% W: matrix with eigenvector columns

[Q,R] = qr(V);

for k = 1:N
    W = A*Q;
    [Q,R] = qr(W);
end
lambda = diag(Q.'*A*Q);
end
