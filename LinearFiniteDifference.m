function [ w ] = LinearFiniteDifference(startp, endp, alpha, beta, N)
%
%LinearFiniteDifference(startp, endp, alpha, beta, N)
%startp, endp: interval
%alpha, beta: y(startp), y(endp), respectively
%N: number of iterations

% y = p(x)y' + q(x)y + r(x)
p = @(x) 0;
q = @(x) 100;
r = @(x) 0;

a = zeros(N+1, 1);
w = a; b = a; c = a; d = a; l = a; u = a; z = a;

h = (endp - startp)/(N + 1);
fprintf('%d', h);
x = startp + h;
a(1) = 2 + q(x)*h^2;
b(1) = (h/2)*p(x) - 1;
d(1) = (1+(h/2)*p(x))*alpha - r(x)*h^2;

for i = 2:(N-1) 
    x = startp + i*h;
    a(i) = 2 + q(x)*h^2;
    b(i) = (h/2)*p(x) - 1;
    c(i) = -(h/2)*p(x) - 1;
    d(i) = -(r(x)*h^2);
end

x = endp - h;
a(N) = 2 + q(x)*h^2;
c(N) = -1 - p(x)*(h/2);
d(N) = (1 - (h/2)*p(x))*beta - r(x)*h^2;

l(1) = a(1);
u(1) = b(1)/a(1);
z(1) = d(1)/l(1);

for i = 2:(N-1)
    l(i) = a(i) - c(i)*u(i-1);
    u(i) = b(i)/l(i);
    z(i) = (d(i) - c(i)*z(i-1))/l(i);
end

l(N) = a(N) - c(N)*u(N-1);
z(N) = (d(N) - c(N)*z(N-1))/l(N);

w(N+1) = beta;
w(N) = z(N);

for i =(N-1):-1:1
    w(i) = z(i) - u(i)*w(i+1);
end

end

