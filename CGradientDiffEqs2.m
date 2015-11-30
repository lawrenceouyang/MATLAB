function [u] = CGradientDiffEqs2(q_xy, r_xy, N, TOL)
%
% CGradientDiffEqs2(q_xy, r_xy, N, TOL)
% This function solves a partial differential equation of the form
% -u_xx - u_yy + q(x,y)u = r(x,y) using the conjugate gradient method.
% q_xy: Function q(x,y) of the equation
% r_xy: Function r(x,y) of the equation
% N: Number of iterations
% TOL: Threshold value for relative residual

h = 1/N;

u = zeros(N-1);
b = zeros(N-1);
q = zeros(N-1);
for i = 1:(N-1)
    for j = 1:(N-1)
        q(i,j) = q_xy(i*h,j*h);
        b(i,j) = h^2*r_xy(i,j);
    end
end

r = zeros(N-1);
gamma = 0;
uxx = 0;
uyy = 0;

for i = 1:(N-1)
    for j = 1:(N-1)
        switch i
            case 1
                uxx = -2*u(i,j) + u(i+1,j);
            case (N-1)
                uxx = u(i-1,j) - 2*u(i,j);
            otherwise
                uxx = u(i-1,j) - 2*u(i,j) + u(i+1,j);
        end
        
        switch j
            case 1
                uyy = -2*u(i,j) + u(i,j+1);
            case (N-1)
                uyy = u(i,j-1) - 2*u(i,j);
            otherwise
                uyy = u(i,j-1) - 2*u(i,j) + u(i,j+1);
        end
        
        r(i,j) = b(i,j) - (-uxx - uyy + h^2*q(i,j)*u(i,j));
    end
end

p = r;
gamma = sum(sum(r.*r));
bbnorm = sum(sum(b.*b));
err = zeros(N-1,1);
a = 0;

w = zeros(N-1);
for k = 1:(N-1)
    
    for i = 1:(N-1)
        for j = 1:(N-1)
        switch i
            case 1
                pxx = -2*p(i,j) + p(i+1,j);
            case (N-1)
                pxx = p(i-1,j) - 2*p(i,j);
            otherwise
                pxx = p(i-1,j) - 2*p(i,j) + p(i+1,j);
        end
        
        switch j
            case 1
                pyy = -2*p(i,j) + p(i,j+1);
            case (N-1)
                pyy = p(i,j-1) - 2*p(i,j);
            otherwise
                pyy = p(i,j-1) - 2*p(i,j) + p(i,j+1);
        end
        
        w(i,j) = -pxx - pyy + h^2*q(i,j)*p(i,j);
            
        end
    end
    
    a = gamma/sum(sum(p.*w));
    u = u + a*p;
    r = r - a*w;
    gamma2 = sum(sum(r.*r));
    p = r + (gamma2/gamma)*p;
    gamma = gamma2;
    err(k) = sqrt(gamma/bbnorm);
    if err(k) < TOL
       fprintf('Solution successfully found after %d iterations.', k);
       err(k+1:N-1) = [];
       loglog(err);
       return;
    end
end
loglog(err);
error('Failed to obtain solution after %d iterations.', sqrt(gamma/bbnorm));
end