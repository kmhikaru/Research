function [L_ord,ord] = comp_L_ord(AAt)
ord = symamd(AAt);
L_ord = chol(AAt(ord,ord),'lower');
