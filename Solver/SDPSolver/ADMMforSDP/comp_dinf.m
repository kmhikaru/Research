function [d_inf,d_err] = comp_dinf(C,S,yA,nob,s)
d_err = norm(C - S - yA);
d_inf = d_err / (1 + sum(abs(C)));
