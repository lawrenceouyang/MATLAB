function [ u, r ] = MatrixRK4( a, b, h, y_a, A )
%
%MatrixRK4( a, b, h, y_a, A )
%Computes the approximate solution of y'(t) = A*y(t) using Runge-Kutta 4
%Inputs:
%a, b: startpoint, endpoint
%h: stepsize
%y_a: initial vector y(a) = y_a
%A: matrix of function

%Outputs:
%u: approximate solution u_N, where N is number of iterations
%r: error calculated using analytic solution y(t) = V*e^(D*t)*V^(-1)*y(a)

[V,D] = eig(A);
A_h = h*A;
N = (b - a)/h;
u = y_a;
err = zeros(N, 1);

for i = 1:N
    K_1 = A_h * u;
    K_2 = A_h * (u + K_1/2);
    K_3 = A_h * (u + K_2/2);
    K_4 = A_h * (u + K_3);
    
    %updated approximation u_i
    u = u + (K_1 + 2*K_2 + 2*K_3 + K_4)/6;
    
    %Exact Solution: Eq.(6) y_i
    y = V * expm(D*(a + i*h))* inv(V) * y_a;
    
    %||y_i - u_i||^2
    err(i) = norm(y - u)^2;
end
r = sqrt((1/N)*sum(err));
fprintf('Error with timestep %d: %d\n', h, r);
end


