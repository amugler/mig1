function dmdt = mdot_fb(t,m,alpha,eta,H0)

dmdt = alpha./(1+(eta./t).^(H0*m));