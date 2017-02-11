function [p_inf,p_err] = comp_pinf(A,X,b)
p_err = norm(A * X - b);
p_inf = p_err / (1 + norm(b));
