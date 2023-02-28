% fitting analysis of deletion control data
% start time based on rightmost position datapoint
% fit activator model (Hill function) in two stages:
% - fit QR data to find IC m0 and decay rate nu
% - fit all data to find H, alpha, eta = K/k
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

% Fit to all data to find H, alpha, eta = K/k
tscale = mean(t);
mscale = mean(m);
Hs = 2:20;
alphas = logspace(0,3,30);
etas = logspace(0,1,30);
tA = linspace(0,4,100);
dt = tA(2)-tA(1);
for i = 1:length(Hs)
    H = Hs(i)
    for j = 1:length(alphas)
        alpha = alphas(j);
        for k = 1:length(etas)
            eta = etas(k);
            mA = exp(-nuopt*tA).*(m0opt + ...
                alpha*dt*cumsum(exp(nuopt*tA)./(1+(eta./tA).^H)));
            for l = 1:length(t)
                xA = ((t(l)-tA)/tscale).^2;
                yA = ((m(l)-mA)/mscale).^2;
                Cl(l) = min(xA + yA);
            end
            CA(i,j,k) = sum(Cl);
        end
    end
end
[CAmin,n] = min(CA(:));
[i,j,k] = ind2sub(size(CA),n);
Hopt = Hs(i)
alphaoptA = alphas(j)
etaopt = etas(k)
mA = exp(-nuopt*tA).*(m0opt + ...
    alphaoptA*dt*cumsum(exp(nuopt*tA)./(1+(etaopt./tA).^Hopt)));

% Plotting
ms = 20; ms2 = 10; ms3 = 13;
lw = 1.5; lw2 = 1;
pu = [.5 0 .5];
gr = .75*[1 1 1];

figure(1); clf
subplot(2,2,1)
hold on
h = plot(tc,mc,'r.',tb,mb,'g.',ta,ma,'b.','markersize',ms);
ht2 = plot(tA,mA,'c-','linewidth',lw);
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
imagesc(Hs,log10(alphas),log10(squeeze(CA(:,:,k)))')
plot(Hopt,log10(alphaoptA),'wo')
xlim([min(Hs) max(Hs)])
ylim([log10(min(alphas)) log10(max(alphas))])
colorbar
xlabel('H')
ylabel('log_{10} \alpha')
title('log_{10} C_A')
set(gca,'ydir','normal','layer','top')
box on

subplot(2,2,3); hold on
imagesc(Hs,log10(etas),log10(squeeze(CA(:,j,:)))')
plot(Hopt,log10(etaopt),'wo')
xlim([min(Hs) max(Hs)])
ylim([min(log10(etas)) max(log10(etas))])
colorbar
xlabel('H')
ylabel('log_{10} K/k')
title('log_{10} C_A')
set(gca,'ydir','normal','layer','top')
box on

subplot(2,2,4); hold on
imagesc(log10(alphas),log10(etas),log10(squeeze(CA(i,:,:)))')
plot(log10(alphaoptA),log10(etaopt),'wo')
xlim([min(log10(alphas)) max(log10(alphas))])
ylim([min(log10(etas)) max(log10(etas))])
colorbar
xlabel('log_{10} \alpha')
ylabel('log_{10} K/k')
title('log_{10} C_A')
set(gca,'ydir','normal','layer','top')
box on

save([dd 'fig23_processed1.mat'])
