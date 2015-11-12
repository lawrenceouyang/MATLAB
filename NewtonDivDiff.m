function [F_x] = NewtonDivDiff(x, f_x)
%Using Algorithm 3.2 on page 126
%x and f_x are both arrays assumed with the same size.
%Due to matlab indices, 1 has to appended to several of
%of the functions.

% Declare array size
X = zeros(length(x));
for i = 0:length(x)-1
    X(i+1,1) = f_x(i+1);
end
for i = 2:length(x)
    for j = 2:i
        X(i,j) = (X(i,j-1)-X(i-1,j-1))/(x(i)-x(i-j+1));
    end
end
F_x = zeros(1,length(x));
for i = 1:length(x)
    F_x(i) = X(i,i);
end

%3.3.7a.
%x = [-0.1, 0, 0.2, 0.3]
%f_x = [5.3, 2, 3.19, 1]
%doing: NewtonDivDiff(x, f_x)
% gives us the result:
%5.3 -33.00, 129.8333, -556.6667
% and the polynomial is:
%p_3(x) = 5.3-33(x+0.1)+129.8333(x+0.1)(x)-5556.6667(x+0.1)(x)(x-0.2)
%
%3.3.7b.
%x = [-0.1, 0, 0.2, 0.3, 0.35]
%f_x = [5.3, 2, 3.19, 1, 0.97260]
%doing: NewtonDivDiff(x, f_x)
% gives us the result:
%5.3 -33.00, 129.8333, -556.6667
% and the polynomial is:
%p_4(x) =
%5.3-33(x+0.1)+129.8333(x+0.1)(x)-5556.6667(x+0.1)(x)(x-0.2)+
%2730.2(x+0.1)(x)(x-0.2)(x-0.3)
