function [X0,S0] = init_X_S(blk,nob,sob)

for j = 1:nob
    if j == 1
      X0 = vec(speye(blk(j)));
    else
      X0 = vertcat(X0,vec(speye(blk(j))));
    end
end
len = length(X0);
S0 = sparse(zeros(len,1));
