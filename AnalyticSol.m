function [ u ] = AnalyticSol( a, b, h, y_a, A )
%
%AnalyticSol( a, b, h, y_a, A )
%Computes the exact solution of y'(t) = Ay(t)
%a, b: startpoint, endpoint
%h: stepsize
%y_a: initial vector y(a) = y_a
%A: matrix of function

%This function redundantly iterates through all steps between a and b of
%step-size h. This was a precursor function used to build RK4 and Adams.
[V,D] = eig(A);
N = (b-a)/h;
u = y_a;

for i = 1:N
    %Exact Solution: Eq.(6) y_i
    u = V*expm(D*(a+i*h))*inv(V)*y_a;
end
end

