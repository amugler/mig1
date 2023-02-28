% fitting analysis of deletion control data
% start time based on rightmost position datapoint
% fit repressor model (Hill function) in two stages:
% - fit QR data to find IC m0 and decay rate nu
% - fit all data to find alpha, beta = (r0/K)^H, gamma = mu*H
% using minimum distance for residuals
clear all

dd = './';
M1 = load([dd 'fig234_data_ctrl.txt']);
for i = 1:5
    ind{i} = find(M1(:,3) == i);
end
pa = M1(ind{1},2);
pb = M1(ind{2},2);
pc = M1(ind{3},2);
p = [pa; pb; pc];
pmax = max(pa);
t = pmax-p;
ta = pmax-pa;
tb = pmax-pb;
tc = pmax-pc;
ma = M1(ind{1},1);
mb = M1(ind{2},1);
mc = M1(ind{3},1);
m = [ma; mb; mc];

% Exponential fit to QR data to find m0 and nu
ya = log(ma);
c = polyfit(ta,ya,1);
nuopt = -c(1)
m0opt = exp(c(2))

% Fit to all data to find alpha, beta = (r0/K)^H, gamma = mu*H
tscale = mean(t);
mscale = mean(m);
gammas = logspace(0,2,30);
alphas = logspace(1,3,30);
betas = logspace(6,9,30);
tR = linspace(0,4,100);
dt = tR(2)-tR(1);
for i = 1:length(gammas)
    gamma = gammas(i)
    for j = 1:length(alphas)
        alpha = alphas(j);
        for k = 1:length(betas)
            beta = betas(k);
            mR = exp(-nuopt*tR).*(m0opt + ...
                alpha*dt*cumsum(exp(nuopt*tR)...
                ./(1+beta*exp(-gamma*tR))));
            for l = 1:length(t)
                xA = ((t(l)-tR)/tscale).^2;
                yA = ((m(l)-mR)/mscale).^2;
                Cl(l) = min(xA + yA);
            end
            CR(i,j,k) = sum(Cl);
        end
    end
end
[CRmin,n] = min(CR(:));
[i,j,k] = ind2sub(size(CR),n);
gammaopt = gammas(i)
alphaoptR = alphas(j)
betaopt = betas(k)
mR = exp(-nuopt*tR).*(m0opt + ...
    alphaoptR*dt*cumsum(exp(nuopt*tR)...
    ./(1+betaopt*exp(-gammaopt*tR))));

% Plotting
ms = 20; ms2 = 10; ms3 = 13;
lw = 1.5; lw2 = 1;
pu = [.5 0 .5];
gr = .75*[1 1 1];
or = [241 162 57]/256;

figure(1); clf
subplot(2,2,1)
hold on
h = plot(tc,mc,'r.',tb,mb,'g.',ta,ma,'b.','markersize',ms);
ht2 = plot(tR,mR,'-','linewidth',lw,'color',or);
xlim([-.2 3.8])
ylim([-1 45])
xlabel('Time, t (AU)')
ylabel('mRNA spots, m')
title('Control')
set(gca,'xdir','reverse')
% legend([h([3 2 1]);hi;hj;hs;ht2],...
%     {'QR','QR.p','QR.pa','least pts','most pts',...
%     '(t*, m*)','degradation'},...
%     'position',[.56 .81 0 0])
box on

subplot(2,2,2); hold on
imagesc(log10(gammas),log10(alphas),log10(squeeze(CR(:,:,k)))')
plot(log10(gammaopt),log10(alphaoptR),'wo')
xlim([log10(min(gammas)) log10(max(gammas))])
ylim([log10(min(alphas)) log10(max(alphas))])
colorbar
xlabel('log_{10} \muH')
ylabel('log_{10} \alpha')
title('log_{10} C_R')
set(gca,'ydir','normal','layer','top')
box on

subplot(2,2,3); hold on
imagesc(log10(gammas),log10(betas),log10(squeeze(CR(:,j,:)))')
plot(log10(gammaopt),log10(betaopt),'wo')
xlim([log10(min(gammas)) log10(max(gammas))])
ylim([min(log10(betas)) max(log10(betas))])
colorbar
xlabel('log_{10} \muH')
ylabel('H log_{10} r_0/K')
title('log_{10} C_R')
set(gca,'ydir','normal','layer','top')
box on

subplot(2,2,4); hold on
imagesc(log10(alphas),log10(betas),log10(squeeze(CR(i,:,:)))')
plot(log10(alphaoptR),log10(betaopt),'wo')
xlim([min(log10(alphas)) max(log10(alphas))])
ylim([min(log10(betas)) max(log10(betas))])
colorbar
xlabel('log_{10} \alpha')
ylabel('H log_{10} r_0/K')
title('log_{10} C_R')
set(gca,'ydir','normal','layer','top')
box on

save([dd 'fig23_processed2.mat'])
