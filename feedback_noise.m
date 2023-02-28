clear all
dd = './';

% plot
figure(1); clf
fs = 12; fs2 = 11;
pu = [148 33 146]/255;

% plot data noise
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
p2c = M2(ind2{3},2);
p3c = M3(ind3{3},2);
t1b = pmax-p1b;
t1c = pmax-p1c;
t1 = [t1b; t1c];
t2b = pmax-p2b;
t2c = pmax-p2c;
t2 = [t2b; t2c];
t3b = pmax-p3b;
t3c = pmax-p3c;
t3 = [t3b; t3c];

% for normalizing by mean QR.pa time in plots
t1Qa = mean(t1c);
t2Qa = mean(t2c);
t3Qa = mean(t3c);
tmax = 1.4;

% get mRNA numbers
m1b = M1(ind1{2},1);
m1c = M1(ind1{3},1);
m1 = [m1b; m1c];
m2b = M2(ind2{2},1);
m2c = M2(ind2{3},1);
m2 = [m2b; m2c];
m3b = M3(ind3{2},1);
m3c = M3(ind3{3},1);
m3 = [m3b; m3c];

% get mean and noise
mthresh = 25;
i1 = find(m1 >= mthresh);
i2 = find(m2 >= mthresh);
i3 = find(m3 >= mthresh);
t1_ = t1(i1);
t2_ = t2(i2);
t3_ = t3(i3);
t1bar = mean(t1_);
t2bar = mean(t2_);
t3bar = mean(t3_);
var_t1 = var(t1_);
var_t2 = var(t2_);
var_t3 = var(t3_);

subplot(3,3,7)
bar([var_t1/t1bar^2 var_t2/t2bar^2 var_t3/t3bar^2],...
    'facecolor',.75*[1 1 1])
ylim([0 .043])
ylabel('CV^2','fontsize',fs)
set(gca,'xtick',[1 2 3],'xticklabels',...
    {'Ctrl','ga80','\DeltaN'},'fontsize',fs2)


% plot data and auto fit
load([dd 'auto_noise.mat'])
lw = 1.5; ms = 12; fs = 12; fs2 = 11;

subplot(3,3,1); hold on
plot(t/t1Qa,a,'r-','linewidth',lw);
set(gca,'xdir','reverse','fontsize',fs2)
xlim([-.05 tmax])
ylim([-1 150])
%xlabel('← Time (AU)','fontsize',fs)
ylabel('Activator number','fontsize',fs)
title('Control','fontsize',fs)

subplot(3,3,2); hold on
plot(t/t2Qa,a,'r-','linewidth',lw);
set(gca,'xdir','reverse','fontsize',fs2)
xlim([-.05 tmax])
ylim([-1 150])
%xlabel('← Time (AU)','fontsize',fs)
%ylabel('Activator number','fontsize',fs)
title('{\it bar-1}(ga80)','fontsize',fs)

subplot(3,3,3); hold on
plot(t/t3Qa,a,'r-','linewidth',lw);
set(gca,'xdir','reverse','fontsize',fs2)
xlim([-.05 tmax])
ylim([-1 150])
%xlabel('← Time (AU)','fontsize',fs)
%ylabel('Activator number','fontsize',fs)
title('\DeltaN-BAR-1','fontsize',fs)

subplot(3,3,4); hold on
plot(t1c/t1Qa,m1c,'m.',t1b/t1Qa,m1b,'g.','markersize',ms);
plot(t/t1Qa,m1,'r-','linewidth',lw);
set(gca,'xdir','reverse','fontsize',fs2)
xlim([-.05 tmax])
ylim([-1 47])
xlabel('← Time (AU)','fontsize',fs)
ylabel('{\it mig-1} number','fontsize',fs)

subplot(3,3,5); hold on
plot(t2c/t2Qa,m2c,'m.',t2b/t2Qa,m2b,'g.','markersize',ms);
plot(t/t2Qa,m2,'r-','linewidth',lw);
set(gca,'xdir','reverse','fontsize',fs2)
xlim([-.05 tmax])
ylim([-1 47])
xlabel('← Time (AU)','fontsize',fs)
%ylabel('{\it mig-1} number','fontsize',fs)

subplot(3,3,6); hold on
plot(t3c/t3Qa,m3c,'m.',t3b/t3Qa,m3b,'g.','markersize',ms);
plot(t/t3Qa,m3,'r-','linewidth',lw);
set(gca,'xdir','reverse','fontsize',fs2)
xlim([-.05 tmax])
ylim([-1 47])
xlabel('← Time (AU)','fontsize',fs)
%ylabel('{\it mig-1} number','fontsize',fs)

% plot auto noise
subplot(3,3,8)
bar([sigmat1/tbar1 sigmat2/tbar2 sigmat3/tbar3].^2,'facecolor','r')
ylim([0 .043])
%ylabel('CV','fontsize',fs)
set(gca,'xtick',[1 2 3],'xticklabels',...
    {'Ctrl','ga80','\DeltaN'},'fontsize',fs2)


% plot both fit
load([dd 'both_noise.mat'])
lw = 1.5;
amin = 10;
amax = 140;

subplot(3,3,1)
plot(t/t1Qa,a1,'-','linewidth',lw,'color',pu);

subplot(3,3,2)
plot([0 1.5],[amin amin],'-','linewidth',lw,'color',pu);

subplot(3,3,3)
plot([0 1.5],[amax amax],'-','linewidth',lw,'color',pu);

subplot(3,3,4)
plot(t/t1Qa,m1,'-','linewidth',lw,'color',pu);

subplot(3,3,5)
plot(t/t2Qa,m2,'-','linewidth',lw,'color',pu);

subplot(3,3,6)
plot(t/t3Qa,m3,'-','linewidth',lw,'color',pu);

% plot both noise
subplot(3,3,9)
bar([sigmat1/tbar1 sigmat2/tbar2 sigmat3/tbar3].^2,'facecolor',pu)
ylim([0 .043])
%ylabel('CV','fontsize',fs)
set(gca,'xtick',[1 2 3],'xticklabels',...
    {'Ctrl','ga80','\DeltaN'},'fontsize',fs2)


% plot cycle fit
subplot(3,3,1)
plot(t/t1Qa,a1,'c--','linewidth',lw);
box on

subplot(3,3,2)
plot([0 1.5],[amin amin],'c--','linewidth',lw);
box on

subplot(3,3,3)
plot([0 1.5],[amax amax],'c--','linewidth',lw);
box on

subplot(3,3,4)
plot(t/t1Qa,m1,'c--','linewidth',lw);
box on

subplot(3,3,5)
m2 = alpha/(1+(K/amin)^H)*t;
plot(t/t2Qa,m2,'c--','linewidth',lw);
box on

subplot(3,3,6)
m3 = alpha/(1+(K/amax)^H)*t;
plot(t/t3Qa,m3,'c--','linewidth',lw);
box on


fd = './';
print(gcf,'-depsc2','-painters','noise.eps')