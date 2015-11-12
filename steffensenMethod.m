function p = steffensenMethod(p_0, TOL, N_0)
%Uses fixed point iteration acceleration method
%to find a solution to p = g(p) given approximation
%p_0 within N_0 iterations and a error of TOL

    function x = g(x)
        x = 2*x;
    end

    i = 1;
    while (i < N_0)
        p_1 = g(p_0);
        p_2 = g(p_1);
        p = p_0 - (p_1 - p_0)^2/(p_2 - 2*p_1 + p_0);
        if (abs(p-p_0)) < TOL)
            return;
        else
            i=i+1;
            p_0 = p;
        end
    end        
    error('The mehtod failed after %d iterations', N_0);
end
