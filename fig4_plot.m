% Activator and mig-1 distributed according to volume fraction
% mig-1 not produced in QR.pp
% taking m(0) = 0
% 
% Activator production constant
% mig-1 production dependent on activator number

clear all
dd = './'; 

% load division time, rho, and model parameters
load([dd 'fig4_processed3.mat']); % activator fit
nu = nuopt;
eta = etaopt; % eta = K/k
H = Hopt;
alpha = alphaoptA;

% example activator and volume parameters for illustration
K = 20;
V0 = 1;
k = K/eta;

% time vectors, incl on either side of division time
t = linspace(0,max(tA),1e3);
ta = linspace(0,T,1e3);
tb = linspace(T,4,1e3);
dta = ta(2)-ta(1);
dtb = tb(2)-tb(1);

% WT size ratio
Z = 4;
rho0 = rhobar;
rhos = linspace(rho0,1,Z);
fs = rhos./(1+rhos);

% QR.p
[V_p,a_p,m_p] = act_dyn(ta,V0,k,0,nu,alpha,K,H,0);
% no division
[V_nd,a_nd,m_nd] = act_dyn(tb-T,V0,k,a_p(end),nu,alpha,K,H,m_p(end));

V_pa = []; a_pa = []; m_pa = [];
for i = 1:Z
    f = fs(i);
    % QR.pa
    [V_pa(:,i),a_pa(:,i),m_pa(:,i)] = act_dyn(tb-T,f*V0,...
        k,f*a_p(end),nu,alpha,K,H,f*m_p(end));    
    % QR.pp
    [V_pp(:,i),a_pp(:,i),m_pp(:,i)] = act_dyn(tb-T,(1-f)*V0,...
        k,(1-f)*a_p(end),nu,0,K,H,(1-f)*m_p(end));
end

% plot properties
lw = 2; lw2 = .5;
fs = 14; fs2 = 12;
gr = .85*[1 1 1];
co = linspace(0,.85,Z);

% plot: V/a/m analysis
F1 = figure(1); clf

% volume
subplot(2,2,1)
hold on
yl = [0 1.3];
plot([T T],yl,'k--','linewidth',lw2)
h(1) = plot(ta,V_p,'g-','linewidth',lw);
h(2) = plot(tb,V_nd,'g--','linewidth',lw);
for i = 1:Z
    h(2+Z+i) = plot(tb,V_pp(:,i),'-',...
        'linewidth',lw,'color',[co(i) 1 1]);
    h(2+i) = plot(tb,V_pa(:,i),'-',...
        'linewidth',lw,'color',[1 co(i) 1]);
end
xlim([-.2 3.8])
ylim(yl)
xlabel('← Time, t (AU)')
ylabel('Volume, V (AU)')
legend(h([1 2 3 3+Z]),...
    {'QR.p','no div.','QR.pa','QR.pp'},'location','se')
set(gca,'xdir','reverse','fontsize',fs)
box on

% activator
subplot(2,2,2)
hold on
yl = [-2 45];
plot([T T],yl,'k--','linewidth',lw2)
plot(ta,a_p,'g-','linewidth',lw);
plot(tb,a_nd,'g--','linewidth',lw);
for i = 1:Z
    plot(tb,a_pp(:,i),'-',...
        'linewidth',lw,'color',[co(i) 1 1]);
    plot(tb,a_pa(:,i),'-',...
        'linewidth',lw,'color',[1 co(i) 1]);
end
xlim([-.2 3.8])
ylim(yl)
xlabel('← Time, t (AU)')
ylabel('Activator number, a')
set(gca,'xdir','reverse','fontsize',fs)
box on

% activator conc.
subplot(2,2,3)
hold on
yl = [-2 45];
plot([T T],yl,'k--','linewidth',lw2)
plot(ta,a_p./V_p,'g-','linewidth',lw);
plot(tb,a_nd./V_nd,'g--','linewidth',lw);
for i = 1:Z
    plot(tb,a_pp(:,i)./V_pp(:,i),'-',...
        'linewidth',lw,'color',[co(i) 1 1]);
    plot(tb,a_pa(:,i)./V_pa(:,i),'-',...
        'linewidth',lw,'color',[1 co(i) 1]);
end
xlim([-.2 3.8])
ylim(yl)
xlabel('← Time, t (AU)')
ylabel('Activator conc., a/V')
set(gca,'xdir','reverse','fontsize',fs)
box on

% mig-1
subplot(2,2,4)
hold on
yl = [-2 45];
plot([T T],yl,'k--','linewidth',lw2)
plot(ta,m_p,'g-','linewidth',lw);
plot(tb,m_nd,'g--','linewidth',lw);
for i = 1:Z
    plot(tb,m_pp(:,i),'-',...
        'linewidth',lw,'color',[co(i) 1 1]);
    plot(tb,m_pa(:,i),'-',...
        'linewidth',lw,'color',[1 co(i) 1]);
end
xlim([-.2 3.8])
ylim(yl)
xlabel('← Time, t (AU)')
ylabel('mig-1 mRNAs, m')
set(gca,'xdir','reverse','fontsize',fs)
box on

% inset: zoom-in on division time
Fi = get(F1,'CurrentAxes');
axes('pos',[.76 .22 .12 .2]);
hold on
plot([T T],yl,'k--','linewidth',lw2)
plot(ta,m_p,'g-','linewidth',lw);
plot(tb,m_nd,'g--','linewidth',lw);
for i = 1:Z
    plot(tb,m_pp(:,i),'-',...
        'linewidth',lw,'color',[co(i) 1 1]);
    plot(tb,m_pa(:,i),'-',...
        'linewidth',lw,'color',[1 co(i) 1]);
end
xlim([2.2 2.45])
ylim([0 7])
set(gca,'xdir','reverse','fontsize',fs)
box on

% plot: to compare with experimental data
figure(2); clf
tV2 = pmax - 1; % time of V2
or = [255 147 0]/255;
gr = .8*[1 1 1];

tmin = 3;
tmid = 3.15;
tmax = 3.3;

% pig-1 null
subplot(2,2,1)
hold on
yl = [-2 47];
plot([tV2 tV2],[9.5 yl(1)],'k:','linewidth',lw)
plot([tV2 tV2],[22 yl(2)],'k:','linewidth',lw)
%plot([tmin tmin],yl,'k-','linewidth',lw)
%plot([tmax tmax],yl,'k-','linewidth',lw)
fill([tmin tmin tmax tmax],[yl yl(2) yl(1)],gr,'linestyle','none')
plot(ta,m_p,'g-','linewidth',lw);
for i = Z:-1:1
    plot(tb,m_pp(:,i),'-',...
        'linewidth',lw,'color',or + co(i)*([1 1 1]-or));
    hpig(i) = plot(tb,m_pa(:,i),'-',...
        'linewidth',lw,'color',[1 co(i) 1]);
    lstr{i} = ['\rho = ' num2str(round(10*rhos(i))/10)];
end
annotation('arrow',[.21 .15],[.71 .65],'linewidth',lw)
text(2.8,16,{'Decreasing','size ratio'},'fontsize',fs)
xlim([0 3.8])
ylim(yl)
xlabel('← Time (AU)')
ylabel('{\it mig-1} mRNAs')
%text(tV2-.1,40,'V2','fontsize',fs)
%legend(hpig,lstr,'location','ne')
set(gca,'xdir','reverse','fontsize',fs)
a = gca; set(a,'box','off','color','none');
b = axes('position',get(a,'position'),'box','on','xtick',[],'ytick',[]);
axes(a);
set(gca,'layer','top')

% m vs rho
rho_ = linspace(.7,3.5,1e2);
f_ = rho_./(1+rho_);
for i = 1:length(f_)
    f = f_(i);
    tbmin = linspace(T,tmin,1e3);
    [V_,a_,mmin(:,i)] = act_dyn(tbmin-T,f*V0,...
        k,f*a_p(end),nu,alpha,K,H,f*m_p(end));
    tbmid = linspace(T,tmid,1e3);
    [V_,a_,mmid(:,i)] = act_dyn(tbmid-T,f*V0,...
        k,f*a_p(end),nu,alpha,K,H,f*m_p(end));
    tbmax = linspace(T,tmax,1e3);
    [V_,a_,mmax(:,i)] = act_dyn(tbmax-T,f*V0,...
        k,f*a_p(end),nu,alpha,K,H,f*m_p(end));
end
subplot(2,2,2)
hold on
fill([rho_ rho_(end:-1:1)],[mmin(end,:) mmax(end,end:-1:1)],...
    gr,'linestyle','none');
plot(rho_,mmid(end,:),'k-','linewidth',lw)
xlim([.8 3.4])
ylim([0 26])
xlabel('QR.pa : QR.pp size ratio')
ylabel('{\it mig-1} mRNAs')
%legend(hpig2,lstr,'location','nw')
set(gca,'fontsize',fs,'xtick',1:.5:3,'ytick',0:5:25)
a = gca; set(a,'box','off','color','none');
b = axes('position',get(a,'position'),'box','on','xtick',[],'ytick',[]);
set(b,'layer','top')
axes(a);

fd = './';
print(gcf,'-depsc',[fd 'fig4EF.eps'])
