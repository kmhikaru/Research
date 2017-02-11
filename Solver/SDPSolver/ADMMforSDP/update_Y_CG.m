function [y] = update_Y_CG(A,b,c,X,y,S,myu,AA,cg_eps)

b1 =  myu * b + A * (c - S - myu * X);

r = b1 - AA*y;
p = r;

for i = 1:length(b) + 1
  if norm(r) < cg_eps
    break;
  end

  alf = (r' * r) / (p' * AA *p);

  y = y + alf * p;

  pre_r = r;
  r = pre_r - alf * AA * p;

  % if norm(r) < cg_eps
  %   break;
  % end

  bt = (r' * r) / (pre_r' * pre_r);
  p = r + bt * p;
end
