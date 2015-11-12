function [ x ] = JacobiIteration(n, A, b, X0, TOL, N)
%
%JacobiIteration(n, A, b, X0, TOL, N)
%n: number of equations/unknowns (max(size(b))
%A: Matrix of the equation Ax = b
%b: b of Ax = b
%X0: initial approximation x0
%TOL: tolerance
%N: maximum number of iterations
x = X0;

for z=1:N
    for i=1:n
        sum = 0;
        for j=1:n
            if i~=j
                sum = sum + A(i,j)*X0(j);
            end
        end
        x(i) = (-sum + b(i))/A(i,i);
    end
    if norm(x - X0) < TOL
        fprintf('Success on iteration #%d',z);
        return;
    end
    X0 = x;
end
fprintf('Maximum number of iterations exceeded.\n');
end