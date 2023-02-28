clear all

dd = './';
load([dd 'auto_lof.mat'])
alphaL = alphaoptA;
etaL = etaopt;
HL = Hopt;
load([dd 'auto_gof.mat'])
alphaG = alphaoptA;
etaG = etaopt;
HG = Hopt;
load([dd 'auto_ctrl.mat'])
alpha = alphaopt;
eta = etaopt;
H0 = H0opt;

k = 30;
KL = k*etaL;
KG = k*etaG;
K = k*eta;
mstar = 25;
T = 4.5;
A = round(k*T + sqrt(10*k*T))

% deterministic (m1: WT, m2: LOF, m3: GOF)
m0 = 0;
t = linspace(0,T+1,1e4);
dt = t(2)-t(1);
[t_out,m1_out] = ode15s(@(t,m) ...
    mdot_auto(t,m,alpha,eta,H0),t,m0);
m1 = m1_out';
m2 = alphaL*dt*cumsum(1./(1+(etaL./t).^HL));
m3 = alphaG*dt*cumsum(1./(1+(etaG./t).^HG));
a = k*t;

figure(1); clf
subplot(2,2,1)
f = 5;
plot([0 T],[A A]/f,'b:',[0 T],mstar*[1 1],'r:',...
    t,a/f,'b-',t,m1,'r-',t,m2,'r--',t,m3,'r-.');
xlim([0 T])
ylim([0 1.1*max([A/f mstar])])
xlabel('t')
ylabel('molecule number')
legend({['A/' num2str(f)],'m*',['a/' num2str(f)],'m (Ctrl)',...
    'm (LOF)','m (GOF)'},'location','best')

% stochastic: numerical
avec = (0:A)';
mvec = (0:mstar-1)';
amat = avec*ones(1,mstar);
mmat = ones(A+1,1)*mvec';
kvec = repmat([k*ones(A,1);0],mstar,1);

% LOF
f = alphaL./(1+(KL./amat).^HL);
fvec = reshape(f,(A+1)*mstar,1);
M0 = -(kvec+fvec);
M_1 = kvec(1:end-1);
M_A1 = fvec(1:end-(A+1));
M = diag(M0) + diag(M_1,-1) + diag(M_A1,-(A+1));
Minv = inv(M);
P0 = zeros((A+1)*mstar,1);
P0(1) = 1;
V = [zeros((A+1)*(mstar-1),1);f(:,end)];
tbar2 = V'*Minv^2*P0;
t2bar = -2*V'*Minv^3*P0;
sigmat2 = sqrt(t2bar-tbar2^2);

% GOF
f = alphaG./(1+(KG./amat).^HG);
fvec = reshape(f,(A+1)*mstar,1);
M0 = -(kvec+fvec);
M_1 = kvec(1:end-1);
M_A1 = fvec(1:end-(A+1));
M = diag(M0) + diag(M_1,-1) + diag(M_A1,-(A+1));
Minv = inv(M);
P0 = zeros((A+1)*mstar,1);
P0(1) = 1;
V = [zeros((A+1)*(mstar-1),1);f(:,end)];
tbar3 = V'*Minv^2*P0;
t2bar = -2*V'*Minv^3*P0;
sigmat3 = sqrt(t2bar-tbar3^2);

% WT
f = alpha./(1+(K./amat).^(H0*mmat));
fvec = reshape(f,(A+1)*mstar,1);
M0 = -(kvec+fvec);
M_1 = kvec(1:end-1);
M_A1 = fvec(1:end-(A+1));
M = diag(M0) + diag(M_1,-1) + diag(M_A1,-(A+1));
Minv = inv(M);
P0 = zeros((A+1)*mstar,1);
P0(1) = 1;
V = [zeros((A+1)*(mstar-1),1);f(:,end)];
tbar1 = V'*Minv^2*P0;
t2bar = -2*V'*Minv^3*P0;
sigmat1 = sqrt(t2bar-tbar1^2);

subplot(2,2,2); hold on
bar([tbar1 tbar2 tbar3])
ylabel('<t>')
box on

subplot(2,2,3); hold on
bar([sigmat1^2 sigmat2^2 sigmat3^2])
ylabel('\sigma_t^2')
box on

subplot(2,2,4); hold on
bar([sigmat1^2/tbar1^2 sigmat2^2/tbar2^2 sigmat3^2/tbar3^2])
ylabel('\sigma_t^2/<t>^2')
box on

save([dd 'auto_noise.mat'])
