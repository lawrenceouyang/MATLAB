function [ A ] = Householder( n, A )
%
% Householder(n, A)
% Obtains symmetric diagonal matrix A(n-1) similar to the symmetric matrix
% A using Householder's method.
% n: dimension
% A: symmetric matrix

for k = 1:(n-2)
    
    q = 0;
    alpha = 0;
    PROD = 0;
    v = zeros(n,1);
    u = zeros(n,1);
    z = zeros(n,1);
    
    for j = (k+1):n
        q = q + A(j,k)^2;
    end
    
    if A(k+1,k) == 0
        alpha = -sqrt(q);
    else
        alpha = -(sqrt(q)*A(k+1,k)/abs(A(k+1,k)));
    end
    
    RSQ = alpha^2 - alpha*A(k+1,k);
    v(k+1) = A(k+1,k) - alpha;
    
    for j = (k+2):n
        v(j) = A(j,k);
    end
    
    for j = k:n
        for i = (k+1):n
            u(j) = u(j) + A(j,i)*v(i);
        end
        u(j) = u(j)/RSQ;
    end
    
    for i = (k+1):n
        PROD = PROD + v(i)*u(i);
    end
    
    for j = k:n
        z(j) = u(j) - (PROD/(2*RSQ))*v(j);
    end
    
    for l = (k+1):n-1
        for j = (l+1):n
            A(j,l) = A(j,l) - v(l)*z(j) - v(j)*z(l);
            A(l,j) = A(j,l);
        end
        A(l,l) = A(l,l) - 2*v(l)*z(l);
    end
    A(n,n) = A(n,n) - 2*v(n)*z(n);
    
    for j = (k+2):n
        A(k,j) = 0;
        A(j,k) = 0;
    end
    A(k+1,k) = A(k+1,k) - v(k+1)*z(k);
    A(k,k+1) = A(k+1,k);
end
end

