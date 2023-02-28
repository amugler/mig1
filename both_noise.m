clear all

dd = './';
load([dd 'both_lof.mat'])
alphaL = alphaLopt;
load([dd 'both_gof.mat'])
alphaG = alphaGopt;
load([dd 'both_ctrl.mat'])
alpha = alphaoptA;
eta = etaopt;
H = Hopt;

beta0 = 25;
K = beta0*eta;
mstar = 25;
J = 10;
beta = 3*beta0;
A = 150;

% deterministic (m1: WT, m2: LOF, m3: GOF)
T = 5;
t = linspace(0,T,1e4);
X0 = [0 0];
[t_out,X_out] = ode15s(@(t,X) ...
    Xdot_cycle(t,X,alpha,K,H,beta0,beta,J),t,X0);
a1 = X_out(:,1);
m1 = X_out(:,2);
m2 = alphaL*t;
m3 = alphaG*t;

figure(1); clf
subplot(2,2,1)
f = 5;
plot([0 T],[A A]/f,'b:',[0 T],mstar*[1 1],'r:',...
    t,a1/f,'b-',t,m1,'r-',t,m2,'r--',t,m3,'r-.');
xlim([0 T])
ylim([0 1.1*max([A/f mstar])])
xlabel('t')
ylabel('molecule number')
legend({['A/' num2str(f)],'m*',['a/' num2str(f) ' (Ctrl)'],...
    'm (Ctrl)','m (LOF)','m (GOF)'},'location','best')

% stochastic: analytic for LOF, GOF
tbar2 = mstar/alphaL;
tbar3 = mstar/alphaG;
sigmat2 = sqrt(mstar)/alphaL;
sigmat3 = sqrt(mstar)/alphaG;

% stochastic: numerical for WT
avec = (0:A)';
mvec = (0:mstar-1)';
amat = avec*ones(1,mstar);
mmat = ones(A+1,1)*mvec';

gmat = beta0 + beta./(1+(J./mmat).^H);
gmat(end,:) = 0;
gVec = reshape(gmat,(A+1)*mstar,1);

fvec = alpha./(1+(K./avec).^H);
fVec = repmat(fvec,mstar,1);

M0 = -(gVec+fVec);
M_1 = gVec(1:end-1);
M_A1 = fVec(1:end-(A+1));
M = diag(M0) + diag(M_1,-1) + diag(M_A1,-(A+1));
Minv = inv(M);

P0 = zeros((A+1)*mstar,1);
P0(1) = 1;
V = [zeros((A+1)*(mstar-1),1);fvec];
tbar1 = V'*Minv^2*P0
t2bar = -2*V'*Minv^3*P0;
sigmat1 = sqrt(t2bar-tbar1^2)

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

save([dd 'both_noise.mat'])
