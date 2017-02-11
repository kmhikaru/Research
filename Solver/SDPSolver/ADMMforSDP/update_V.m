function [V_next,yA] = update_V(A,C,X,Y,myu)
%normal update
yA = (Y' * A)';
V_next = C - yA - myu * X;
