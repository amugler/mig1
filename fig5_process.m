% variance analysis of bar-1 data
% Brown-Forsythe ANOVA test
% start time based on rightmost QR.p position of all strains
% 1. t: time = const - position
% 2. s: time rescaled by mean QR.pa time
% 3. r: time rescaled by mean time of all points above mthresh

clear all

dd = './';
M1 = load([dd 'fig5_data_ctrl.txt']);
M2 = load([dd 'fig5_data_gof.txt']);
M3 = load([dd 'fig5_data_lof.txt']);
for i = 1:5
    ind1{i} = find(M1(:,3) == i);
    ind2{i} = find(M2(:,3) == i);
    ind3{i} = find(M3(:,3) == i);
end

% find rightmost QR.p position and define time t
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

% define rescaled time s
s1b = t1b/mean(t1c);
s1c = t1c/mean(t1c);
s1 = [s1b; s1c];
s2b = t2b/mean(t2c);
s2c = t2c/mean(t2c);
s2 = [s2b; s2c];
s3b = t3b/mean(t3c);
s3c = t3c/mean(t3c);
s3 = [s3b; s3c];

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

% variance and ANOVA test
mmin = 10;
mmax = 25;
mthresh = mmin:mmax;
for k = 1:length(mthresh)
    i1 = find(m1 >= mthresh(k));
    i2 = find(m2 >= mthresh(k));
    i3 = find(m3 >= mthresh(k));
    c1 = ones(length(i1),1);
    c2 = 2*ones(length(i2),1);
    c3 = 3*ones(length(i3),1);

    % t
    t1_ = t1(i1);
    t2_ = t2(i2);
    t3_ = t3(i3);
    var_t1(k) = var(t1_);
    var_t2(k) = var(t2_);
    var_t3(k) = var(t3_);
    pCGt(k) = vartestn([t1_;t2_],[c1;c2],...
        'testtype','brownforsythe','display','off');
    pCLt(k) = vartestn([t1_;t3_],[c1;c3],...
        'testtype','brownforsythe','display','off');
    p3t(k) = vartestn([t1_;t2_;t3_],[c1;c2;c3],...
        'testtype','brownforsythe','display','off');

    % s
    s1_ = s1(i1);
    s2_ = s2(i2);
    s3_ = s3(i3);
    var_s1(k) = var(s1_);
    var_s2(k) = var(s2_);
    var_s3(k) = var(s3_);
    pCGs(k) = vartestn([s1_;s2_],[c1;c2],...
        'testtype','brownforsythe','display','off');
    pCLs(k) = vartestn([s1_;s3_],[c1;c3],...
        'testtype','brownforsythe','display','off');
    p3s(k) = vartestn([s1_;s2_;s3_],[c1;c2;c3],...
        'testtype','brownforsythe','display','off');

    % r
    r1_ = t1_/mean(t1_);
    r2_ = t2_/mean(t2_);
    r3_ = t3_/mean(t3_);
    var_r1(k) = var(r1_);
    var_r2(k) = var(r2_);
    var_r3(k) = var(r3_);
    pCGr(k) = vartestn([r1_;r2_],[c1;c2],...
        'testtype','brownforsythe','display','off');
    pCLr(k) = vartestn([r1_;r3_],[c1;c3],...
        'testtype','brownforsythe','display','off');
    p3r(k) = vartestn([r1_;r2_;r3_],[c1;c2;c3],...
        'testtype','brownforsythe','display','off');
end

% plot
figure(1); clf
lw = 1; ms = 5; fs = 11;

% t
% plot variances vs threshold
subplot(3,3,1)
plot(mthresh,var_t1,'k.-',...
    mthresh,var_t3,'m.-',...
    mthresh,var_t2,'c.-','linewidth',lw,'markersize',ms)
xlabel('m threshold')
ylabel('Variance')
title('Position')
%legend({'Control','LOF','GOF'},'location','best')
set(gca,'fontsize',fs)

% plot p-values
subplot(3,3,2)
h = semilogy([mmin mmax],.05*[1 1],'k:',...
    mthresh,pCGt,'b.-',mthresh,2*pCGt,'b--',...
    mthresh,pCLt,'r.-',mthresh,2*pCLt,'r--',...
    mthresh,p3t,'g.-','linewidth',lw,'markersize',ms);
xlabel('m threshold')
ylabel('p value')
% legend(h(2:end),{'Control-GOF','Bonferoni-corrected (x2)',...
%     'Control-LOF','Bonferoni-corrected (x2)',...
%     'All 3'},'location','se')
set(gca,'fontsize',fs)

% plot variances for a particular threshold
subplot(3,3,3)
mstar = 10;
k = find(mthresh==mstar);
bar([var_t1(k) var_t3(k) var_t2(k)],'facecolor',.75*[1 1 1])
set(gca,'xticklabels',{'Ctrl','LOF','GOF'})
ylabel('Variance')
title(['m >= ' num2str(mstar)])
set(gca,'fontsize',fs)

% s
% plot variances vs threshold
subplot(3,3,4)
plot(mthresh,var_s1,'k.-',...
    mthresh,var_s3,'m.-',...
    mthresh,var_s2,'c.-','linewidth',lw,'markersize',ms)
xlabel('m threshold')
ylabel('Variance')
title('Rescaled by <QR.pa>')
%legend({'Control','LOF','GOF'},'location','best')
set(gca,'fontsize',fs)

% plot p-values
subplot(3,3,5)
h = semilogy([mmin mmax],.05*[1 1],'k:',...
    mthresh,pCGs,'b.-',mthresh,2*pCGs,'b--',...
    mthresh,pCLs,'r.-',mthresh,2*pCLs,'r--',...
    mthresh,p3s,'g.-','linewidth',lw,'markersize',ms);
xlabel('m threshold')
ylabel('p value')
% legend(h(2:end),{'Control-GOF','Bonferoni-corrected (x2)',...
%     'Control-LOF','Bonferoni-corrected (x2)',...
%     'All 3'},'location','se')
set(gca,'fontsize',fs)

% plot variances for a particular threshold
subplot(3,3,6)
mstar = 25;
k = find(mthresh==mstar);
bar([var_s1(k) var_s3(k) var_s2(k)],'facecolor',.75*[1 1 1])
set(gca,'xticklabels',{'Ctrl','LOF','GOF'})
ylabel('Variance')
title(['m >= ' num2str(mstar)])
set(gca,'fontsize',fs)

% r
% plot variances vs threshold
subplot(3,3,7)
plot(mthresh,var_r1,'k.-',...
    mthresh,var_r3,'m.-',...
    mthresh,var_r2,'c.-','linewidth',lw,'markersize',ms)
xlabel('m threshold')
ylabel('Variance')
title('Rescaled by <t>')
%legend({'Control','LOF','GOF'},'location','best')
set(gca,'fontsize',fs)

% plot p-values
subplot(3,3,8)
h = semilogy([mmin mmax],.05*[1 1],'k:',...
    mthresh,pCGr,'b.-',mthresh,2*pCGr,'b--',...
    mthresh,pCLr,'r.-',mthresh,2*pCLr,'r--',...
    mthresh,p3r,'g.-','linewidth',lw,'markersize',ms);
xlabel('m threshold')
ylabel('p value')
% legend(h(2:end),{'Control-GOF','Bonferoni-corrected (x2)',...
%     'Control-LOF','Bonferoni-corrected (x2)',...
%     'All 3'},'location','se')
set(gca,'fontsize',fs)

% plot variances for a particular threshold
subplot(3,3,9)
mstar = 25;
k = find(mthresh==mstar);
bar([var_r1(k) var_r3(k) var_r2(k)],'facecolor',.75*[1 1 1])
set(gca,'xticklabels',{'Ctrl','LOF','GOF'})
ylabel('Variance')
title(['m >= ' num2str(mstar)])
set(gca,'fontsize',fs)

pCGmax = max(2*pCGr)
pCGmin = min(2*pCGr)
pCGmean = mean(2*pCGr)
pCLmax = max(2*pCLr)
pCLmin = min(2*pCLr)
pCLmean = mean(2*pCLr)
