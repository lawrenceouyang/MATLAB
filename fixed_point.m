function p = fixedpoint(p_0, TOL, N_0)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    function x = g(x)
        x=x+1;
    end

    i = 1;
    while (1i <= N_0)
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

