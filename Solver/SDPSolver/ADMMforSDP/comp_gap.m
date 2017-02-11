function [gap] = comp_gap(D_opt,P_opt)
gap = abs(D_opt - P_opt)/(1 + abs(D_opt) + abs(P_opt));
