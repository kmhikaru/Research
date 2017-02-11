%Author : Hikaru
%this method is Alternating direction method of multipliers for Linear Programming
%this program solves LP such as below.
%====================LP====================
% maximize      c^T x
% subject to    Ax <= b
%               x >= 0
%==========================================
% File format
% --------m: size of b, n: size of c
% m n
% b1 b2 b3 ... bm
% c1 c2 c3 ... cn
% a11 a12 ...a1n
% ....
% an1 an2 ...amn

function optval = admmForLP(fileName)
optval = '#';
fileID = fopen(fileName, 'r');
formatSpec = '%f';
MN = fscanf(fileID, '%d', 2);
M = MN(1);
N = MN(2);

b = fscanf(fileID, formatSpec, M);
c = -fscanf(fileID, formatSpec, N);
A = fscanf(fileID, formatSpec, [N,M]);
A = transpose(A);

B = horzcat(A,eye(M));
Bt = transpose(B);
BBt = B*Bt;
BBt_inv = inv(BBt);
p = vertcat(c,zeros(M,1));
%==== initialize z and s ====
z = ones(N+M,1);
s = zeros(N+M,1);
w = ones(M,1);
pre_z = z;





cnt = 2000000;

tol = 1e-7;
check_gap = 1;

myu = 10;
rho = 1.6;

%======= sub parameter to change main parameter myu =====
eta1 = 10;
eta2 = 10;
itr_pinf = 0;
itr_dinf = 0;
h = 1000;
gmm = 0.7;
myu_min = 1e-5;
myu_max = 10000;

for i = 1:cnt
  %======= updating variavles =============
  w = -BBt_inv*(B*(myu*z + s - p) - myu * b);
  %w = cg(z,s,w,B,BBt,b,p,myu);
  %normyL = B*z - b + B*(Bt*w + s - p)/myu;
  v = p - Bt*w - myu*z;
  s = diag(v)*(v>0);
  z = (s - v)/myu;
  z = (1-rho)*pre_z + rho*z; %momentum
  pre_z = z;
  %========================================


  p_opt = p'*z;
  d_opt = b'*w;
  gap = abs(p_opt - d_opt)/(1+abs(p_opt)+abs(d_opt));
  pinf = norm(B*z - b)/(1 + norm(b));
  dinf = norm(p - s - B'*w);


  %check condition
  if max([gap,pinf,dinf]) < tol
    optval = -p_opt;
    % fprintf('%e \n',optval);
    break;
  end
  if i > 10000 & (abs(p_opt/d_opt) > 10 | abs(d_opt/p_opt) > 10)
    % fprintf('#\n');
    optval = '#';
    break;
  end


  if i == cnt & max([gap,pinf,dinf]) < 0.1
    % optval = -(p_opt+d_opt)/2;
    optval = -p_opt;
  end

  %see the　progress
  % if max([gap,pinf,dinf]) < check_gap
    %  fprintf('Iter: %d primal opt: %e dual opt: %e pinf: %e dinf: %e gap: %e myu: %e \n',i,p_opt,d_opt,pinf,dinf,gap,myu);
  %   check_gap = check_gap*0.1;
  %   if gap < tol
  %     fprintf('%e \n',p_opt);
  %     break;
  %   end
  % end

  %change parameters

  if i == 100000
    h = 30000;
    gmm = 0.9999;
    tol = tol * 10;
  end
  if i == 300000
    tol = tol *10;
  end
  if pinf/dinf <= eta1
    itr_pinf = itr_pinf +1;
    itr_dinf = 0;
    if itr_pinf >= h
      myu = max([gmm*myu,myu_min]);
      itr_pinf = 0;
    end
  end
  if pinf/dinf > eta2
    itr_dinf = itr_dinf + 1;
    itr_pinf = 0;
    if itr_dinf >= h
      myu = min([myu/gmm,myu_max]);
      itr_dinf = 0;
    end
  end

end
