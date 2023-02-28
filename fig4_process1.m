% finding the division time from the control data

clear all

dd = './';
M1 = load([dd 'fig234_data_ctrl.txt']);
for i = 1:5
    ind{i} = find(M1(:,3) == i);
end
pa = M1(ind{1},2);
pb = M1(ind{2},2);
pc = M1(ind{3},2);
pmax = max(pa);
tb = pmax-pb;
tc = pmax-pc;

tb_sort = sort(tb,'descend');
Tb = length(tb);
tc_sort = sort(tc);
Tc = length(tc);

% intersection time is 5th (and 6th) time in sorted tc
T = tc_sort(5)

figure(1); clf; hold on
lw = 2;
stairs(tb_sort,1:Tb,'g','linewidth',lw)
stairs(tc_sort,1:Tc,'m','linewidth',lw)
plot([T T],[0 12],'k--','linewidth',lw)
xlim([-.2 3.8])
ylim([0 12])
xlabel('Time (AU)')
ylabel('Cumulative number of data points')
set(gca,'xdir','reverse')
box on

save([dd 'fig4_processed1.mat'])
