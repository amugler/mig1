% fitting analysis of LOF data: constant activation with rate alphaL
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
p2c = M2(ind2{3},2);
t2b = pmax-p2b;
t2c = pmax-p2c;
t2 = [t2b; t2c];

% get mRNA numbers
m2b = M2(ind2{2},1);
m2c = M2(ind2{3},1);
m2 = [m2b; m2c];

% Fit to QR.p, QR.pa data to find alphaL
tscale = mean(t2);
mscale = mean(m2);
alphaLs = logspace(0,2,100);
tA = linspace(0,5,1e3);
for j = 1:length(alphaLs)
    alphaL = alphaLs(j);
    mA = alphaL*tA;
    for l = 1:length(t2)
        xA = ((t2(l)-tA)/tscale).^2;
        yA = ((m2(l)-mA)/mscale).^2;
        Cl(l) = min(xA + yA);
    end
    CA(j) = sum(Cl);
end
[CAmin,j] = min(CA);
alphaLopt = alphaLs(j)
mA = alphaLopt*tA;

% Plotting
ms = 20; ms2 = 10; ms3 = 13;
lw = 1.5; lw2 = 1;
pu = [.5 0 .5];
gr = .75*[1 1 1];

figure(1); clf
subplot(2,2,1)
hold on
h = plot(t2c,m2c,'m.',t2b,m2b,'g.','markersize',ms);
ht2 = plot(tA,mA,'c-','linewidth',lw);
xlim([-.2 4.5])
ylim([-1 45])
xlabel('Time, t (AU)')
ylabel('mRNA spots, m')
title('LOF')
set(gca,'xdir','reverse')
box on

subplot(2,2,2)
loglog(alphaLs,CA)
xlabel('\alpha_L')
ylabel('C_A')

save([dd 'both_lof.mat'])
