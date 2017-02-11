function [x,y,s,data] = admm_r(A,b,c,K,params,parCoLO)

%Author : tetteketei
%this method is Alternating direction method of multipliers for Semidefinite Programing
%this program solves LP such as below.
%====================LP====================
% maximize      c^T x
% subject to    Ax = b
%               x \in K ... K is SDP cone
%==========================================


% params struct consists of the following fields
% here set to default settings
max_iters = 5000; % maximum num iterations for admm
eps = 1e-3;        % quitting tolerances
cg_eps = 1e-13;    % quitting tolerances in CG method
use_indirect = 2;  % 0 ... direct / 1 ... use CG method / 2 ... use sparse Cholesky Factorization
mu = 5;            % step size


% params for updating %%%%%%%%%
h = 100;
eta1 = 10;
eta2 = 10;
gmm = 0.99;
mu_min = 0.0001;
mu_max = 1000000;

it_p_inf = 0;
it_d_inf = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

res_name = 'res.txt';

if nargin >= 5
  if isfield(params,'max_iters'); max_iters=params.max_iters;end;
  if isfield(params,'eps'); eps=params.eps;end;
  if isfield(params,'cg_eps'); cg_eps=params.cg_eps;end;
  if isfield(params,'mu'); mu=params.mu;end;
  if isfield(params,'use_indirect'); use_indirect=params.use_indirect;end;
  if isfield(params,'h'); h=params.h;end;
  if isfield(params,'eta1'); eta1=params.eta1;end;
  if isfield(params,'eta2'); eta2=params.eta2;end;
  if isfield(params,'gmm'); gmm=params.gmm;end;
  if isfield(params,'mu_min'); mu_min=params.mu_min;end;
  if isfield(params,'mu_max'); mu_max=params.mu_max;end;
  if isfield(params,'res_name'); res_name=params.res_name;end;
end

if nargin == 6
  parCoLO.SDPsolver = [];
  J.f = size(A,1);
  if parCoLO.EQorLMI == 1
    [x,y,infoCoLO,cliqueDomain,cliqueRange,LOP] = sparseCoLO(A,b,c,K,J,parCoLO);
    A = LOP.A;
    b = LOP.b;
    c = LOP.c;
    K = LOP.K;
  elseif parCoLO.EQorLMI == 2
    [x,y,infoCoLO,cliqueDomain,cliqueRange,LOP] = sparseCoLO(A,b,c,K,J,parCoLO);
    [A,b,c,K] = lmi2eq(LOP.A,LOP.b,LOP.c,LOP.J);
  end
end

AAt  = A*A';
% [A,b] = scaling(A,b);
% cond(AAt)

blk = K.s;
nob = length(blk);
sob = length(b);

for j = 1:nob
  len(j) = blk(j)^2;
end

% pre allocation
l = sum(len);
X_next = zeros(l,1);
S_next = zeros(l,1);
V_next = zeros(l,1);
Y_next = ones(sob,1);

scl = 0.1;
ak=1;

fprintf('\nSettings ........ eps: %e gmm: %e h: %d eta1: %d eta2: %d \n',eps,gmm,h,eta1,eta2);

% initializing
fprintf('Initializing ');
[X_next,S_next] = init_X_S(blk,nob,sob);
fprintf('.......=>Finish \n');


if use_indirect == 0
  fprintf('Computing AAt and its inverse ');
  tic
  % AAt = A*A';
  AAt_inv = sparse(inv(AAt));
  fprintf('.......=>Finish \n');
  toc
elseif use_indirect == 1
  % AAt  = A*A';
elseif use_indirect == 2
  fprintf('Sparse Cholesky Factorization ');
  tic
  [L_ord,ord] = comp_L_ord(A*A');
  Lt_ord = L_ord';
  fprintf('.......=>Finish \n');
  toc
end

tic

for i = 1:max_iters
  if use_indirect == 0 % direct method
    Y_next = update_Y_Inv(A,b,c,AAt_inv,X_next,S_next,mu); %analytical solution
  elseif use_indirect == 1 % use CG method
    Y_next = update_Y_CG(A,b,c,X_next,Y_next,S_next,mu,AAt,cg_eps); %CG method
  elseif use_indirect == 2 % use sparse Cholesky Factorization
    Y_next = update_Y_LL(L_ord,Lt_ord,ord,A,b,c,X_next,S_next,mu)';
  end

  data.err_y(i) = norm(AAt * Y_next + mu * (A * X_next - b) + A * (S_next - c));
  normyL = norm(A * X_next - b + (A * (A' * Y_next + S_next - c )) / mu);

  [V_next,yA] = update_V(A,c,X_next,Y_next,mu);

  X_pre = X_next;

  [S_next,X_next] = update_S_X(X_pre,V_next,mu,nob,blk);

  %=========== momentum============
  ak1 = (1 + sqrt(1 + 4 * ak * ak)) / 2;
  rho = (ak1 + ak - 1) / ak1;
  [S_next,X_next] = momentum_S_X(X_pre,S_next,V_next,mu,rho);
  ak = ak1;
  clear V_next;
  d_opt = b' * Y_next + 0;
  p_opt = c' * X_next + 0;

  [p_inf,p_err] = comp_pinf(A,X_next,b);
  [d_inf,d_err] = comp_dinf(c,S_next,yA,nob,blk);

  gap = comp_gap(d_opt,p_opt);
  dlt = max([p_inf,d_inf,gap]);

  data.p_inf(i) = p_inf;
  data.p_err(i) = p_err;
  data.d_inf(i) = d_inf;
  data.d_err(i) = d_err;
  data.cx(i) = p_opt;
  data.by(i) = d_opt;
  data.gap(i) = gap;
  data.mu(i) = mu;
  data.normyL(i) = normyL;

  if dlt <= eps
    fprintf('========== dlt < eps:%e ========== \n',eps);
    print_res(i,p_opt,d_opt,gap,mu,p_inf,d_inf,normyL);
    break
  end

  if dlt <= scl
    fprintf('========== dlt < %e ==========\n',scl);
    print_res(i,p_opt,d_opt,gap,mu,p_inf,d_inf,normyL);
    scl = scl * 0.1;
    toc
  end


  [mu,it_p_inf,it_d_inf] = update_mu(p_inf,d_inf,eta1,eta2,it_p_inf,it_d_inf,h,gmm,mu,mu_min,mu_max);

end

x = X_next;
y = Y_next;
s = S_next;
toc
write2Data(res_name,data.p_inf,data.d_inf,data.cx,data.by,data.gap,data.mu,data.normyL);
