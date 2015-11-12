function p = fixedpoint(p_0, TOL, N_0)
%Fixed Point Iteration
%   This function iterates N_0 times the formula defined
%   by g(x) initially using p_0 and then checks if the 
%   error is less than the tolerance TOL. 
%
%   Problem 2.2.8:
%   g(x) = 2^(-x) in interval [1/3,1] and is a decreasing
%   function. Thus max is g(1/3)= 0.7937 and min is  g(1)= 1/2,
%   which satisifies (i) of Theorem 2.3.
%   In addition, abs(g'(x)) = abs(-ln2(2^-x)), is a decreasing
%   function, with value ranging from [0.3465, 0.5501] for x in 
%   [1/3,1] which satisfies (ii), thus showing there is a unique
%   fixed point on the interval [1/3,1]. 
%
%   The result of fixedpoint(1,10^(-4),20) gives us an answer
%   of 0.6412. Modifying the function to return the variable 'i'
%   instead, it is found that at the 12th iteration the solution
%   is found.
%
%   Using corollary 2.5, then abs(p_n - p) <= k^n(2/3) <= 10^(-4),
%   where k = 0.5501. Solving for n gives n rounded up to be >= 15.

    function x = g(x)   %there was no specification for fnc declaration
        x = 501.0625 - 201.0625* exp(-.4*x);     %so I left it nested in to be changed manually.
    end

    i = 1;
    while (i <= N_0)
        p = g(p_0);
        if (abs(p-p_0) < TOL)
            return;
        else
            i=i+1;
            p_0 = p;
        end
    end
    error('The method failed after %d iterations', N_0);
end

