function [u] = CGradientDiffEqs(q_xy, r_xy, N, TOL)
%
% CGradientDiffEqs(q_xy, r_xy, N, TOL)
% This function solves a partial differential equation of the form
% -u_xx - u_yy + q(x,y)u = r(x,y) using the conjugate gradient method.
% q_xy: Function q(x,y) of the equation
% r_xy: Function r(x,y) of the equation
% N: Number of iterations
% TOL: Threshold value for relative residual

% Initialize the step size
h = 1/N;

% Initialize u, w, q, b, r, p, matrix for the iteration
% Initialize error vector and gamma
u = zeros(N-1); 
w = zeros(N-1); q = w; b = q; r = b;
err = zeros(N-1,1);

for i = 1:(N-1)
    for j = 1:(N-1)
        q(i,j) = q_xy(i*h,j*h);
        b(i,j) = (h^2)*r_xy(i,j);
    end
end

for i = 1:(N-1)
    for j = 1:(N-1)
        % Consider boundary conditions
        switch i
            case 1
                u_xx = -2*u(i,j) + u(i+1,j);
            case (N-1)
                u_xx = u(i-1,j) - 2*u(i,j);
            otherwise
                u_xx = u(i-1,j) - 2*u(i,j) + u(i+1,j);
        end
        
        switch j
            case 1
                u_yy = -2*u(i,j) + u(i,j+1);
            case (N-1)
                u_yy = u(i,j-1) - 2*u(i,j);
            otherwise
                u_yy = u(i,j-1) - 2*u(i,j) + u(i,j+1);
        end
        
        r(i,j) = b(i,j) - (-u_xx -u_yy + (h^2)*q(i,j)*u(i,j));
    end
end

p = r;
gamma = sum(sum(r.*r));
ip_bb = sum(sum(b.*b));
a = 0;
ip_pw = 0;

%Begin the algorithm
for k = 1:(N-1)
    
    for i = 1:(N-1)
        for j = 1:(N-1)
            % Consider boundary conditions
            switch i
                case 1
                    p_xx = -2*p(i,j) + p(i+1,j);
                case (N-1)
                    p_xx = p(i-1,j) - 2*p(i,j);
                otherwise
                    p_xx = p(i-1,j) - 2*p(i,j) + p(i+1,j);
            end
            switch j
                case 1
                    p_yy = -2*p(i,j) + p(i,j+1);
                case (N-1)
                    p_yy = p(i,j-1) - 2*p(i,j);
                otherwise
                    p_yy = p(i,j-1) - 2*p(i,j) + p(i,j+1);
            end

            w(i,j) = -p_xx - p_yy + h^2*q(i,j)*p(i,j);
        end
    end
    
    ip_pw = sum(sum(p.*w));
    a = gamma/ip_pw;
    u = u + a*p;
    r = r - a*w;
    newGamma = 0;
    for i = 1:(N-1)
        for j = 1:(N-1)
            newGamma = newGamma + r(i,j)^2;
        end
    end
    p = r + (newGamma/gamma)*p;
    gamma = newGamma;
    err(k) = sqrt(gamma/ip_bb);
    if err(k) < TOL
        fprintf('Solution successfully found after %d iterations.', k);
        err(k+1:N-1) = [];
        loglog(err);
        return;
    end
end
loglog(err);
error('Failed to obtain solution after %d iterations.', sqrt(gamma/ip_bb));
end

