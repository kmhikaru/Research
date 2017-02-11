function [w] = cg(z,s,w,B,BBt,b,p,myu)

b1 = myu * b + B * (p - s - myu * z);

r = b1 - BBt * w;
q = r;

cg_eps = 0.000001;
for i = 1:length(b) + 1

  alf = (r'*r)/(q'*BBt*q);

  w = w + alf * q;

  pre_r = r;
  r = pre_r - alf * BBt * q;

  if norm(r) < cg_eps
    break;
  end

  bt = (r'*r)/(pre_r' * pre_r);
  q = r + bt*q;
end
