% fitting analysis of Ctrl data: standard activator model
%   (as simplified limit of cycle + autoregulation model)
% start time based on rightmost position datapoint of all strains
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

% Fit to QR.p, QR.pa data to find H, alpha, eta = K/k
tscale = mean(t1);
mscale = mean(m1);
Hs = 2:20;
alphas = logspace(2,4,30);
etas = logspace(-1,1,30);
tA = linspace(0,5,1e3);
dt = tA(2)-tA(1);
for i = 1:length(Hs)
    H = Hs(i)
    for j = 1:length(alphas)
        alpha = alphas(j);
        for k = 1:length(etas)
            eta = etas(k);
            mA = alpha*dt*cumsum(1./(1+(eta./tA).^H));
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
Hopt = Hs(i)
alphaoptA = alphas(j)
etaopt = etas(k)
mA = alphaoptA*dt*cumsum(1./(1+(etaopt./tA).^Hopt));

% Plotting
ms = 20; ms2 = 10; ms3 = 13;
lw = 1.5; lw2 = 1;
pu = [.5 0 .5];
gr = .75*[1 1 1];

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

save([dd 'both_ctrl.mat'])
