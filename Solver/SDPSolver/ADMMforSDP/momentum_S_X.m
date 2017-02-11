function [S_next,X_next] = momentum_S_X(X_pre,S_next,V_next,myu,rho)
X_next = (1 - rho) * X_pre + rho * (S_next - V_next) / myu; %accelerate
