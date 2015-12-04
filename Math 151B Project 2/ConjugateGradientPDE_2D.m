function [u, k] = ConjugateGradientPDE_2D(q_xy, r_xy, N, TOL, Iteration)
%
% CGradientDiffEqs(q_xy, r_xy, N, TOL)
% This function solves a partial differential equation of the form
% -u_xx - u_yy + q(x,y)u = r(x,y) using the conjugate gradient method.
% Our initial guess u is defined as an N-1 x N-1 zero matrix.
% The function returns the matrix u, which is the solution, and the integer
% k, which is the number of iterations to achieve the solution.
% q_xy: Function q(x,y) of the equation
% r_xy: Function r(x,y) of the equation
% N: Size of partition of domain
% TOL: Threshold value for relative residual
% Iteration: Maximum number of iterations until program terminates

% Initialize the step size
h = 1/N;

% Initialize u, w, q, b, r, p, matrix for the iteration
% Initialize relative residual vector and gamma
u = zeros(N-1); 
w = zeros(N-1); q = w; b = q; r = b;
err = zeros(Iteration,1);

for i = 1:(N-1)
    for j = 1:(N-1)
        q(i,j) = q_xy(i*h,j*h);
        b(i,j) = (h^2)*r_xy(i,j);
    end
end

%Initialize r = b - Au_0, where A is the operation v = -u_xx - u_yy + h^2q(i,j)u(i,j)
for i = 1:(N-1)
    for j = 1:(N-1)
        % Consider boundary conditions
        switch i
            case 1
                u_xx = -2*u(i,j) + u(i+1,j); % u(0,j) = 0
            case (N-1)
                u_xx = u(i-1,j) - 2*u(i,j); % u(N,j) = 0
            otherwise
                u_xx = u(i-1,j) - 2*u(i,j) + u(i+1,j);
        end
        
        switch j
            case 1
                u_yy = -2*u(i,j) + u(i,j+1); % u(i,0) = 0
            case (N-1)
                u_yy = u(i,j-1) - 2*u(i,j); % u(i,N) = 0
            otherwise
                u_yy = u(i,j-1) - 2*u(i,j) + u(i,j+1);
        end
        
        r(i,j) = b(i,j) - (-u_xx -u_yy + (h^2)*q(i,j)*u(i,j));
    end
end

p = r;
Gam = sum(sum(r.*r));
ip_bb = sum(sum(b.*b));


%Begin the algorithm
for k = 1:Iteration; 
    %Solve w = Ap_k, where A is the operation w = -p_xx - p_yy +
    %h^2q(i,j)p(i,j)
    for i = 1:(N-1)
        for j = 1:(N-1)
            % Consider boundary conditions
            switch i
                case 1
                    p_xx = -2*p(i,j) + p(i+1,j); % p(0,j) = 0
                case (N-1)
                    p_xx = p(i-1,j) - 2*p(i,j); % p(N,j) = 0
                otherwise
                    p_xx = p(i-1,j) - 2*p(i,j) + p(i+1,j);
            end
            % Consider boundary conditions
            switch j
                case 1
                    p_yy = -2*p(i,j) + p(i,j+1); % p(i,0) = 0
                case (N-1)
                    p_yy = p(i,j-1) - 2*p(i,j); % p(i,N) = 0
                otherwise
                    p_yy = p(i,j-1) - 2*p(i,j) + p(i,j+1);
            end

            w(i,j) = -p_xx - p_yy + h^2*q(i,j)*p(i,j);
        end
    end
    % solve for a_k
    ip_pw = sum(sum(p.*w));
    a = Gam/ip_pw;
    
    % update u,r,p and gamma
    u = u + a*p;
    r = r - a*w;
    newGamma = sum(sum(r.*r));
    p = r + (newGamma/Gam)*p;
    Gam = newGamma;
    
    % check if relative residual is less than threshold
    err(k) = sqrt(Gam/ip_bb);
    if err(k) < TOL
        fprintf('Solution successfully found after %d iterations.\n', k);
        err(k+1:N-1) = [];
        loglog(err);
        str = sprintf('Log Scale Residual to Iterations, N = %d', N);
        title(str);
        xlabel('Number of Iterations');
        ylabel('Relative Residual');
        return;
    end
end
error('Failed to obtain solution after %d iterations.', k);
end

