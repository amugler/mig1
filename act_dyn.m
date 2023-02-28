function [V,a,m] = act_dyn(t,V0,k,a0,nu,alpha,K,H,m0)

dt = t(2)-t(1);
V = V0*ones(1,length(t));
a = a0 + k*t;
m = exp(-nu*t).*(m0 + alpha*dt*cumsum(exp(nu*t)./(1+(K./a).^H)));
