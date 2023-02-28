function dXdt = Xdot_cycle(t,X,alpha,K,H,beta0,beta,J)

a = X(1);
m = X(2);
dadt = beta0 + beta./(1+(J./m).^H);
dmdt = alpha./(1+(K./a).^H);
dXdt = [dadt;dmdt];