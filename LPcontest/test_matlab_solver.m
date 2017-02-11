%usage .... matlab -nodisplay -r test_matlab_solver\(\'admmForLP\',\'02\'\);
function [] = test_matlab_solver(solver_name,pro_num_str) %solver_name function must have an only argument.

  pro_dir = strcat('./problems/',pro_num_str,'/');
  solver = str2func(solver_name);
  cd(pro_dir);
  problemList = dir(['*.in']);
  len = length(problemList);
  cnt = 0;
  fprintf('============start============\n')
  wrong_probs = [];

  len = 10;

  for i = 1:len
    fprintf('%d\n',i);
    problem = problemList(i).name;
    %tic
    opt1 = solver(problem); %user solver returns optimal value:
    %toc
    %tic
    opt2 = solveLP(problem);        %linprog returns optimal value:
    %toc
    if opt1 == '#'
      fprintf('your opt: #\n');
    else
      fprintf('your opt: %e\n',opt1);
    end
    if opt2 == '#'
      fprintf('answer  : #\n');
    else
      fprintf('answer  : %e\n',opt2);
    end
    % abs(opt1 - opt2)/(1+abs(opt2))
    if abs(opt1 - opt2)/(1+abs(opt2)) < 0.0001 || opt1 == opt2
    % if abs(opt1 - opt2)/(1+abs(opt2)) < 0.0001 || opt1 == opt2
    % if abs(opt1 - opt2) < 0.0001 || opt1 == opt2
      fprintf('Result-----------------> OK\n');
      cnt = cnt + 1;
    else
      wrong_probs = vertcat(wrong_probs,[i]);
      fprintf('NG\n');
    end
  end
  if cnt == len
    fprintf('congratulations!!');
  else
    fprintf('your score %d / %d',cnt,len);
  end
  wrong_probs
exit
