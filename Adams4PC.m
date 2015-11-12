function [ u_N, r] = Adams4PC( a, b, h, y_a, A )
%
%Adams4PC( a, b, h, y_a, A )
%Computes the approximate solution of y'(t) = A*y(t) using Adams 4th Order
%Predictor Corrector, and using RK4 for the first 3 approximations.
%Inputs:
%a, b: startpoint, endpoint
%h: stepsize
%y_a: initial vector y(a) = y_a
%A_t: matrix of function

%Outputs:
%u: approximate solution u_N, where N is number of iterations
%err: error calculated using analytic solution y(t) = V*e^(D*t)*V'*y(a)

[V,D] = eig(A);
N = (b - a)/h;
u = zeros(size(y_a, 1),N);
u(:, 1) = y_a;
err = zeros(1, N);

%Pre-Multiplied Matrices:
A_59 = 59*A;
A_55 = 55*A;
A_37 = 37*A;
A_19 = 19*A;
A_9 = 9*A;
A_5 = 5*A;

%RK4 for u_1, u_2 and u_3
for i = 1:3
    K_1 = h * A * u(:, i);
    K_2 = h * A * (u(:, i) + K_1/2);
    K_3 = h * A * (u(:, i) + K_2/2);
    K_4 = h * A * (u(:, i) + K_3);
    
    u(:, i+1) = u(:, i) + (K_1 + 2*K_2 + 2*K_3 + K_4)/6;
    
    %Exact Solution: Eq.(6) y_i
    y = V * expm(D*(a + i*h)) * inv(V) * y_a;
    
    %||y_i - u_i||^2
    err(i) = norm(y - u(:, i+1))^2;
end

for i = 4:N
    %predictor:
    u_p = u(:, i) + (h/24)*(A_55*u(:, i) - A_59*u(:, i-1) + A_37*u(:, i-2) - A_9*u(:, i-3));
    
    %corrector:
    u(:, i+1) = u(:, i) + (h/24)*(A_9*u_p + A_19*u(:, i) - A_5*u(:, i-1) + A*u(:, i-2));
    
    %Exact Solution: Eq.(6) y_i
    y = V * expm(D*(a + i*h)) * inv(V) * y_a;
    
    %||y_i - u_i||^2
    err(i) = norm(y - u(:, i+1))^2;
end

r = sqrt((1/N)*sum(err));
fprintf('Error with timestep %d: %d\n', h, r);

u_N = u(:, N+1);

end


