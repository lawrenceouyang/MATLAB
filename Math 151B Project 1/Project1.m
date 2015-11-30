clear;

%Assuming data file is in same directory.
load('proj1_151b_data.mat');

%Declare vectors for holding results
RK4A1 = zeros(11,1);
RK4A2 = zeros(11,1);
RK4A3 = zeros(11,1);
ADPCA1 = zeros(11,1);
ADPCA2 = zeros(11,1);
ADPCA3 = zeros(11,1);

c = zeros(11,1);

%Apply all methods for step sizes h = 1, 1/2, 1/4, 1/8, ..., 1/2^10
for i = 1:11
    c(i) = 1/(2^(i-1));
    [~, RK4A1(i)] = MatrixRK4(0,100,1/(2^(i-1)),y0,A1);
    [~, RK4A2(i)] = MatrixRK4(0,100,1/(2^(i-1)),y0,A2);
    [~, RK4A3(i)] = MatrixRK4(0,100,1/(2^(i-1)),y0,A3);
    [~, ADPCA1(i)] = Adams4PC(0,100,1/(2^(i-1)),y0,A1); 
    [~, ADPCA2(i)] = Adams4PC(0,100,1/(2^(i-1)),y0,A1); 
    [~, ADPCA3(i)] = Adams4PC(0,100,1/(2^(i-1)),y0,A1); 
end

%plot the graphs for the reasonable values
subplot(2,3,1)
loglog(c(3:11), RK4A1(3:11))
title('RK4 - A1')

subplot(2,3,2)
loglog(c(3:11), RK4A2(3:11))
title('RK4 - A2')

subplot(2,3,3)
loglog(c(3:11), RK4A3(3:11))
title('RK4 - A3')

subplot(2,3,4)
loglog(c(4:11), ADPCA1(4:11))
title('APC - A1')

subplot(2,3,5)
loglog(c(4:11), ADPCA2(4:11))
title('APC - A2')

subplot(2,3,6)
loglog(c(4:11), ADPCA3(4:11))
title('APC - A3')

%Fit the non-infinite and reasonable values into a linear function.
polyfit(log(c(3:11)),log(RK4A1(3:11)),1)
polyfit(log(c(3:11)),log(RK4A2(3:11)),1)
polyfit(log(c(3:11)),log(RK4A3(3:11)),1)
polyfit(log(c(3:11)),log(ADPCA1(4:11)),1)
polyfit(log(c(3:11)),log(ADPCA2(4:11)),1)
polyfit(log(c(3:11)),log(ADPCA3(4:11)),1)
