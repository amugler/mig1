% fitting analysis of GOF data: constant activation with rate alphaG
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
p3c = M3(ind3{3},2);
t3b = pmax-p3b;
t3c = pmax-p3c;
t3 = [t3b; t3c];

% get mRNA numbers
m3b = M3(ind3{2},1);
m3c = M3(ind3{3},1);
m3 = [m3b; m3c];

% Fit to QR.p, QR.pa data to find alphaG
tscale = mean(t3);
mscale = mean(m3);
alphaGs = logspace(0,2,100);
tA = linspace(0,5,1e3);
for j = 1:length(alphaGs)
    alphaG = alphaGs(j);
    mA = alphaG*tA;
    for l = 1:length(t3)
        xA = ((t3(l)-tA)/tscale).^2;
        yA = ((m3(l)-mA)/mscale).^2;
        Cl(l) = min(xA + yA);
    end
    CA(j) = sum(Cl);
end
[CAmin,j] = min(CA);
alphaGopt = alphaGs(j)
mA = alphaGopt*tA;

% Plotting
ms = 20; ms2 = 10; ms3 = 13;
lw = 1.5; lw2 = 1;
pu = [.5 0 .5];
gr = .75*[1 1 1];

figure(1); clf
subplot(2,2,1)
hold on
h = plot(t3c,m3c,'m.',t3b,m3b,'g.','markersize',ms);
ht2 = plot(tA,mA,'c-','linewidth',lw);
xlim([-.2 4.5])
ylim([-1 45])
xlabel('Time, t (AU)')
ylabel('mRNA spots, m')
title('GOF')
set(gca,'xdir','reverse')
box on

subplot(2,2,2)
loglog(alphaGs,CA)
xlabel('\alpha_G')
ylabel('C_A')

save([dd 'both_gof.mat'])
