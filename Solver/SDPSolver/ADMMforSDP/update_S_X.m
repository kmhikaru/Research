function [S_next,X_next] = update_S_X(X_pre,V_next,myu,nob,s)
%function [S_next,X_next,eigS] = update_SandX7(V_next,eigS,myu,nob,s)
%this is the fastest 12/25
X_next = zeros([sum(s.^2),1]);
S_next = zeros([sum(s.^2),1]);

st = 1;
en = 0;

for j = 1:nob
  en = en + s(j) * s(j);
  matVj = mat(V_next(st:en));
  non_sprs = full(matVj);
  [v,D] = eig(non_sprs);
  d = diag(D);

  if length(find(d>0)) < s(j)/2
    S_next(st:en) = vec(v(:,find(d>0)) * diag(d(find(d>0))) * v(:,find(d>0))');
    X_next(st:en) = vec(S_next(st:en) - V_next(st:en)) / myu;
  else
    X_next(st:en) = vec(v(:,find(d<0)) * diag(d(find(d<0))) * v(:,find(d<0))') / (-myu);
    S_next(st:en) = vec(V_next(st:en) + myu * X_next(st:en));
  end

  st = en + 1;

end
