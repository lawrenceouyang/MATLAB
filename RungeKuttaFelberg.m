function [t, w, h] = RungeKuttaFelberg( f_t_y, a, b, init, TOL, hmax, hmin)
%
%RungeKuttaFelberg( f_t_y, a, b, init, TOL, hmax, hmin)
%f_t_y: Function being solved
%a, b: Start and end point
%init: Initial value
%TOL: Tolerance
%hmax: Maximum step size
%hmin: Minimum step size

t = a;
w = init;
h = hmax;
FLAG = 1;
fprintf('%f, %f\n', t, w)

while (FLAG == 1)
    K_1 = h*f_t_y(t,w);
    K_2 = h*f_t_y(t + h/4, w + K_1/4);
    K_3 = h*f_t_y(t + 3*h/8, w + 3*K_1/32 + 9*K_2/32);
    K_4 = h*f_t_y(t + 12*h/13, w + (1932*K_1 - 7200*K_2 + 7296*K_3)/2197);
    K_5 = h*f_t_y(t + h, w + 439*K_1/216 - 8*K_2 + 3680*K_3/513 - 845*K_4/4104);
    K_6 = h*f_t_y(t + h/2, w - 8*K_1/27 + 2*K_2 - 3544*K_3/2565 + 1859*K_4/4104 - 11*K_5/40);
    
    R = abs(K_1/360 - 128*K_3/4275 - 2197*K_4/75240 + K_5/50 + 2*K_6/55)/h;
    if (R <= TOL)
        t = t + h;  
        w = w + 25*K_1/216 + 1408*K_3/2565 + 2197*K_4/4104 - K_5/5;
        fprintf('%f, %f, %f\n', t, w, h);
        %Use the print command below for comparison to actual solution
        %fprintf('%f, %f, %f, %f\n', t, w, h, %actual solution);
    end
    
    delta = 0.84*(TOL/R)^(1/4);
    if (delta <= 0.1)
        h = 0.1*h;
    elseif (delta >= 4)
        h = 4*h;
    else
        h = 8*h;
    end
    
    if (h > hmax)
        h = hmax;
    end
    
    if (t >= b)
        FLAG = 0;
    elseif (t + h > b)
        h = b - t;
    elseif (h < hmin)
        FLAG = 0;
        fprintf('minimum h(%d) exceeded\n', h);
        error('Procedure completed unsuccessfully.');
    end
end
end

%5.5: 1(c):
%fprime = @(t, w) 1 + t/w
%RungeKuttaFelberg(fprime, 1, 2, 2, 10^-4, 0.25, 0.05)
%Resulting Output:
% i |    t_i       w_i       h_i       y_i
%---+----------+---------+---------+---------
% 0 | 1.000000, 2.000000
% 1 | 1.250000, 2.778930, 0.250000, 2.778929
% 2 | 1.500000, 3.608198, 0.250000, 3.608198
% 3 | 1.750000, 4.479329, 0.250000, 4.479328
% 4 | 2.000000, 5.386296, 0.250000, 5.386294

%5.5: 1(d):
%fprime = @(t, w) cos(2*t) + sin(3*t)
%RungeKuttaFelberg(fprime, 0, 1, 1, 10^-4, 0.25, 0.05)
%Resulting Output:
% i |    t_i       w_i       h_i       y_i
%---+----------+---------+---------+---------
% 0 | 0.000000, 1.000000
% 1 | 0.250000, 1.329148, 0.250000, 1.329150
% 2 | 0.500000, 1.730486, 0.250000, 1.730490
% 3 | 0.750000, 2.041467, 0.250000, 2.041472
% 4 | 1.000000, 2.117975, 0.250000, 2.117980

