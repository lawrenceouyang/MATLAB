function [ x ] = GaussSeidel(n, A, b, X0, TOL, N)
%
%GaussSeidel(n, A, b, X0, TOL, N)
%n: number of equations/unknowns (max(size(b))
%A: Matrix of the equation Ax = b
%b: b of Ax = b
%X0: initial approximation x0
%TOL: tolerance
%N: maximum number of iterations
x = X0;

for z=1:N
    for i=1:n
        sumax = 0;
        sumaxo = 0;
        for j=1:(i-1)
           sumax = sumax + A(i,j)*x(j);
        end
        for j =(i+1):n
            sumaxo = sumaxo + A(i,j)*X0(j);
        end
        x(i) = (- sumax - sumaxo + b(i))/A(i,i);   
    end
    if norm(x - X0) < TOL
        fprintf('Success on iteration #%d',z);
        return;
    end
    X0 = x;
end
fprintf('Maximum number of iterations exceeded.\n');
end