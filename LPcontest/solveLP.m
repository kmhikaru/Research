function optval = solveLP(filename)

fileID = fopen(filename, 'r');
formatSpec = '%f';

MN = fscanf(fileID, '%d', 2);
M = MN(1);
N = MN(2);

b = fscanf(fileID, formatSpec, M);
c = fscanf(fileID, formatSpec, N);
A = fscanf(fileID, formatSpec, [N,M]);
A = transpose(A);

% options = optimset('Display','iter','TolFun',1e-10);

options = optimset('Display','off');

[~, optval, exitflag] = linprog(-c, A, b, [], [], zeros(N,1), [], [], options);
  if exitflag == -3
      optval = '#';
      % fprintf('lin:  %s \n',optval);
  else
      % linprog is minimizing function.
      optval = -optval;
      % fprintf('lin:  %e\n',optval);
  end
end
