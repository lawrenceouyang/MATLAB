function [H_x] = HermiteInterpolate(x, f_x, derf_x)

Q = zeros(length(x));
Z = zeros(1,length(x));

for i = 1:length(x)
    Z(2*i-1) = x(i);
    Z(2*i) = x(i);
    Q(2*i-1,1) = f_x(i);
    Q(2*i,1) = f_x(i);
    Q(2*i,2) = derf_x(i);
    if (i ~= 1)
        Q(2*i-1,2) = (Q(2*i-1,1)-Q(2*i-2,1))/(Z(2*i-1)-Z(2*i-2));
    end
end

for i = 3:2*length(x)+1
    for j = 3:i
        Q(i,j) = (Q(i-1,j-2)-Q(i-2,j-2))/(Z(i-1)-Z(i-j+1));
    end
end
H_x = zeros(1,length(x));
for i = 1:2*length(x)+1
    H_x(i) = Q(i,i);
end

%3.4.1: Using HermiteInterpolate(x, f_x, derf_x), where
%       x =      [-0.5, -0.25, 0]
%       f_x =    [-0.02475, 0.3349375, 1.101]
%       derf_x = [0.751, 2.189, 4.002]
%       gives the result:
%       [-0.0248, 0.7510, 2.7510, 1, 0, 0]
%       which gives the polynomial:
%       H(x) = -0.0248 + 0.7510(x+0.5) + 2.7510(x+0.5)^2 +
%       (x+0.5)^2(x+0.25)



