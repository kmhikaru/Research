function [] = print_res(itr,P_opt,D_opt,gap,myu,pinf,dinf,normyL)
fprintf('Iter %d C*X: %e D_opt: %e gap: %e ',itr,P_opt,D_opt,gap );
fprintf('myu: %e pinf: %e dinf: %e normyL: %e \n',myu,pinf,dinf,normyL);
