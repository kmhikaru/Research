function [Y_next] =  update_Y_LL(L_ord,Lt_ord,ord,A,b,c,X,S,myu)
% (L_ord,ord,A',c,b,X_next,S_next,myu)
b1 =  myu * b + A * (c - S - myu * X);
b_ord = b1(ord);
z_ord = L_ord \ b_ord;
Y_next(ord) = Lt_ord \ z_ord;
