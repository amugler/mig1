% Plotting
% comparing fit vs zero initial mig-1
% comparing fit vs less regulator
clear all

dd = './';
load([dd 'fig23_processed1.mat']); % activator fit
load([dd 'fig23_processed2.mat']); % repressor fit

% regulators
% sample values of bare parameters
KA = 20;
KR = 5;
HR = Hopt;
% infer parameters for regulator dynamics
k = KA/etaopt;
r0 = KR*betaopt^(1/HR);
mu = gammaopt/HR;
% regulator dynamics
a = k*tA;
r = r0*exp(-mu*tR);

% zero initial mig-1
mA1 = exp(-nuopt*tA).*(0 + ...
    alphaoptA*dt*cumsum(exp(nuopt*tA)./(1+(etaopt./tA).^Hopt)));
mR1 = exp(-nuopt*tR).*(0 + ...
    alphaoptR*dt*cumsum(exp(nuopt*tR)...
    ./(1+betaopt*exp(-gammaopt*tR))));

% less activator: decrease k
ka = .85*k; etaa = KA/ka;
mA2a = exp(-nuopt*tA).*(m0opt + ...
    alphaoptA*dt*cumsum(exp(nuopt*tA)./(1+(etaa./tA).^Hopt)));
kb = .7*k; etab = KA/kb;
mA2b = exp(-nuopt*tA).*(m0opt + ...
    alphaoptA*dt*cumsum(exp(nuopt*tA)./(1+(etab./tA).^Hopt)));
mA2c = m0opt*exp(-nuopt*tA);

% less repressor: decrease r0
r0a = .7*r0; betaa = (r0a/KR)^HR;
mR2a = exp(-nuopt*tR).*(m0opt + ...
    alphaoptR*dt*cumsum(exp(nuopt*tR)...
    ./(1+betaa*exp(-gammaopt*tR))));
r0b = .4*r0; betab = (r0b/KR)^HR;
mR2b = exp(-nuopt*tR).*(m0opt + ...
    alphaoptR*dt*cumsum(exp(nuopt*tR)...
    ./(1+betab*exp(-gammaopt*tR))));
mR2c = exp(-nuopt*tR).*(m0opt + alphaoptR*(exp(nuopt*tR)-1));

% plot properties
lw = 3; lw2 = 2;
fs = 14; fs2 = 14;
dc = .5*[0 1 1];
dr = .5*[1 0 0];
gr = .85*[1 1 1];

% Fig 1: comparing fit vs zero initial mig-1
figure(1); clf
Y = 20;

subplot(Y,4,[1:4:Y*2+1 2:4:Y*2+2])
hold on
hA1(1) = plot(tA,a,'-','linewidth',lw,'color',gr);
hA1(2) = plot(tA,mA,'c-','linewidth',lw);
hA1(3) = plot(tA,mA1,'--','linewidth',lw,'color',dc);
xlim([-.2 3.8])
ylim([-2 45])
xlabel('                                                      ← Time (AU)')
ylabel('Molecule number')
set(gca,'xdir','reverse','ytick',0:10:40,'fontsize',fs)
legend(hA1,{'Activator','{\it mig-1}, fit',...
    '{\it mig-1}, {\it m}(0) = 0'},...
    'fontsize',fs2,'location','ne')
a = gca; set(a,'box','off','color','none');
b = axes('position',get(a,'position'),'box','on','xtick',[],'ytick',[]);
axes(a);

subplot(Y,4,[3:4:Y*2+3 4:4:Y*2+4])
hold on
hR1(1) = plot(tR,r,'-','linewidth',lw,'color',gr);
hR1(2) = plot(tR,mR,'r-','linewidth',lw);
hR1(3) = plot(tR,mR1,'--','linewidth',lw,'color',dr);
xlim([-.2 3.8])
ylim([-2 45])
%xlabel('← Time, t (AU)')
%ylabel('mig-1 mRNAs, m')
set(gca,'xdir','reverse','ytick',0:10:40,'fontsize',fs)
legend(hR1,{'Repressor','{\it mig-1}, fit',...
    '{\it mig-1}, {\it m}(0) = 0'},...
    'fontsize',fs2,'location','ne')
a = gca; set(a,'box','off','color','none');
b = axes('position',get(a,'position'),'box','on','xtick',[],'ytick',[]);
axes(a);

fd = './';
print(gcf,'-depsc',[fd 'fig2D.eps'])

% Fig 2: comparing fit vs less regulator
figure(2); clf
subplot(Y,4,[1:4:Y*2+1 2:4:Y*2+2])
hold on
%hA2(1) = plot(tA,tA*0,'-','linewidth',lw,'color',gr);
plot(tA,mA,'c-','linewidth',lw);
plot(tA,mA2a,'-','linewidth',lw,'color',5/6*[0 1 1]);
plot(tA,mA2b,'-','linewidth',lw,'color',2/3*[0 1 1]);
plot(tA,mA2c,'-','linewidth',lw,'color',dc);
% P1 = [2.5 20]; P2 = [3.5 2]; dP = P2-P1;
% quiver(P1(1),P1(2),dP(1),dP(2),0);
annotation('arrow',[.21 .15],[.61 .53],'linewidth',lw2)
text(2.7,13,{'Less','activator'},'fontsize',fs)
xlim([-.2 3.8])
ylim([-2 45])
xlabel('                                                      ← Time (AU)')
ylabel('{\it mig-1} mRNAs')
set(gca,'xdir','reverse','ytick',0:10:40,'fontsize',fs)
% legend(hA2(2:3),{'mig-1 (fit to Control data)',...
%     'mig-1 (with no activator)'},...
%     'fontsize',fs2,'position',[.25 .1 eps eps])
a = gca; set(a,'box','off','color','none');
b = axes('position',get(a,'position'),'box','on','xtick',[],'ytick',[]);
axes(a);


subplot(Y,4,[3:4:Y*2+3 4:4:Y*2+4])
hold on
%hR2(1) = plot(tR,tR*0,'-','linewidth',lw,'color',gr);
plot(tR,mR,'r-','linewidth',lw);
plot(tR,mR2a,'-','linewidth',lw,'color',5/6*[1 0 0]);
plot(tR,mR2b,'-','linewidth',lw,'color',2/3*[1 0 0]);
plot(tR,mR2c,'-','linewidth',lw,'color',dr);
annotation('arrow',[.625 .87],[.61 .83],'linewidth',lw2)
text(3.65,10,{'Less','repressor'},'fontsize',fs)
xlim([-.2 3.8])
ylim([-2 45])
% xlabel('← Time, t (AU)')
% ylabel('mig-1 mRNAs, m')
set(gca,'xdir','reverse','ytick',0:10:40,'fontsize',fs)
% legend(hR2(2:3),{'mig-1 (fit to Control data)',...
%     'mig-1 (with no repressor)'},...
%     'fontsize',fs2,'position',[.75 .1 eps eps])
a = gca; set(a,'box','off','color','none');
b = axes('position',get(a,'position'),'box','on','xtick',[],'ytick',[]);
axes(a);

fd = './';
print(gcf,'-depsc',[fd 'fig3A.eps'])
