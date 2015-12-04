% project2main.m
% This script will run the correct operations as to the desire of the specs.

% Call ConjugateGradientPDE_2D for N = 32, 64, 128, 256, 512, which will
% return iteration numbers and draw plots.

% Initialize our q(x,y) and r(x,y):
q_xy = @(x,y) exp(x+y);
r_xy = @(x,y) 1;

% Consider a maximum number of iterations to run and a threshold for our
% residual:
Iterations = 2000;
TOL = 10^-6;

IterationTot = zeros(5, 2);

%Compute the solutions and the graphs:
for i = 1:5
    N = 2^(i+4);
    IterationTot(i,1) = N;
    [temp, IterationTot(i,2)] = ConjugateGradientPDE_2D(q_xy, r_xy, N, TOL, Iterations);
    str = sprintf('logplotN-%d.jpg',N);
    saveas(gcf,str);
    close;
end

%Compute graph of iterations to N
loglog(IterationTot(1:5,1), IterationTot(1:5,2));
title('Log Graph of N and k');
xlabel('N');
ylabel('k');
saveas(gcf,'logplotN-k.jpg');
close;