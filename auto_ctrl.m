% fitting analysis of Control data: H that increases with m
% start time based on rightmost position datapoint
% using minimum distance for residuals
clear all

dd = './';
M1 = load([dd 'fig5_data_ctrl.txt']);
M2 = load([dd 'fig5_data_lof.txt']);
M3 = load([dd 'fig5_data_gof.txt']);
for i = 1:5
    ind1{i} = find(M1(:,3) == i);
    ind2{i} = find(M2(:,3) == i);
    ind3{i} = find(M3(:,3) == i);
end

% find rightmost QR.p position and define time
p1b = M1(ind1{2},2);
p2b = M2(ind2{2},2);
p3b = M3(ind3{2},2);
pmax = max([p1b; p2b; p3b]);
p1c = M1(ind1{3},2);
t1b = pmax-p1b;
t1c = pmax-p1c;
t1 = [t1b; t1c];

% get mRNA numbers
m1b = M1(ind1{2},1);
m1c = M1(ind1{3},1);
m1 = [m1b; m1c];

% Plotting
ms = 20; ms2 = 10; ms3 = 13;
lw = 1.5; lw2 = 1;
pu = [.5 0 .5];
gr = .75*[1 1 1];

% Fit to QR.p, QR.pa data to find alpha, eta, H0
m0 = 0;
tscale = mean(t1);
mscale = mean(m1);
Z = 20;
alphas = logspace(2,3,Z);
etas = logspace(0,1,Z);
H0s = linspace(1,10,Z);
tA = linspace(0,4,1e3);
for i = 1:length(alphas)
    alpha = alphas(i)
    for j = 1:length(etas)
        eta = etas(j);
        for k = 1:length(H0s)
            H0 = H0s(k);
            [tA_out,mA_out] = ode15s(@(t,m) ...
                mdot_auto(t,m,alpha,eta,H0),tA,m0);
            mA = mA_out';
            for l = 1:length(t1)
                xA = ((t1(l)-tA)/tscale).^2;
                yA = ((m1(l)-mA)/mscale).^2;
                Cl(l) = min(xA + yA);
            end
            CA(i,j,k) = sum(Cl);
        end
    end
end
[CAmin,n] = min(CA(:));
CAmin
[i,j,k] = ind2sub(size(CA),n);
alphaopt = alphas(i)
etaopt = etas(j)
H0opt = H0s(k)
[tA_out,mA] = ode15s(@(t,m) mdot_auto(t,m,alphaopt,etaopt,H0opt),tA,m0);

% plot fit
figure(1); clf
subplot(2,2,1)
hold on
h = plot(t1c,m1c,'m.',t1b,m1b,'g.','markersize',ms);
ht2 = plot(tA,mA,'c-','linewidth',lw);
xlim([-.2 4.5])
ylim([-1 45])
xlabel('Time, t (AU)')
ylabel('mRNA spots, m')
title('Ctrl')
set(gca,'xdir','reverse')
box on

subplot(2,2,2); hold on
imagesc(log10(alphas),log10(etas),log10(squeeze(CA(:,:,k)))')
plot(log10(alphaopt),log10(etaopt),'wo')
xlim([min(log10(alphas)) max(log10(alphas))])
ylim([min(log10(etas)) max(log10(etas))])
colorbar
xlabel('log_{10} \alpha')
ylabel('log_{10} \eta')
title('log_{10} C_A')
set(gca,'ydir','normal','layer','top')
box on

subplot(2,2,3); hold on
imagesc(log10(alphas),H0s,log10(squeeze(CA(:,j,:)))')
plot(log10(alphaopt),H0opt,'wo')
xlim([min(log10(alphas)) max(log10(alphas))])
ylim([min(H0s) max(H0s)])
colorbar
xlabel('log_{10} \alpha')
ylabel('H_0')
title('log_{10} C_A')
set(gca,'ydir','normal','layer','top')
box on

subplot(2,2,4); hold on
imagesc(log10(etas),H0s,log10(squeeze(CA(i,:,:)))')
plot(log10(etaopt),H0opt,'wo')
xlim([min(log10(etas)) max(log10(etas))])
ylim([min(H0s) max(H0s)])
colorbar
xlabel('log_{10} \eta')
ylabel('H_0')
title('log_{10} C_A')
set(gca,'ydir','normal','layer','top')
box on

save([dd 'auto_ctrl.mat'])