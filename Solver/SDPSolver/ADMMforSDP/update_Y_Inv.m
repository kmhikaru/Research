function[y] = update_Y_Inv(A,b,c,AAt_inv,X,S,myu)
y = - AAt_inv * (myu * (A * X - b) + A * (S - c)); 
