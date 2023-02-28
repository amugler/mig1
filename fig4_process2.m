% finding the mean size ratio

clear all

dd = './';
M1 = load([dd 'fig4_data_ratio.txt']);
rhobar = mean(M1)
drho = std(M1)

save([dd 'fig4_processed2.mat'])
