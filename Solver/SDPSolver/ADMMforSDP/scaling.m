function [At,b] = scaling(At,b)

AAt = At*At';
h1 = At.*At;
sv = sqrt(sum(h1')');

for i = 1:length(AAt)
  At(i,:) = At(i,:)./sv(i,1);
end

b = b./sv;
