clear;
close all;
addpath('sim_data/');



load('sim_data_start_2021_07_09_15:48_end_2021_07_09_15:54_sim_total_corr_pure_ICA.mat');%P=2
data_america(1,:) = america_mean_trial_time(:,1,2,end);   
data_GLICA(1,:) = GLICA_mean_trial_time(:,1,2,end);   
data_QICA(1,:) = QICA_mean_trial_time(:,1,2,end);   
data_sa4ica(1,:) = sa4ica_mean_trial_time(:,1,2,end);   


load('sim_data_start_2021_07_09_22:13_end_2021_07_10_12:25_sim_total_corr_pure_ICA.mat');%P=3
data_america(2,:) = america_mean_trial_time(:,1,2,end);   
data_GLICA(2,:) = GLICA_mean_trial_time(:,1,2,end);   
data_QICA(2,:) = QICA_mean_trial_time(:,1,2,end);   
data_sa4ica(2,:) = sa4ica_mean_trial_time(:,1,2,end);   

load('sim_data_start_2021_07_10_12:29_end_2021_07_12_15:57_sim_total_corr_pure_ICA.mat');%P=5
data_america(3,:) = america_mean_trial_time(:,1,end,end);   
data_GLICA(3,:) = GLICA_mean_trial_time(:,1,end,end);   
data_QICA(3,:) = QICA_mean_trial_time(:,1,end,end);   
data_sa4ica(3,:) = sa4ica_mean_trial_time(:,1,end,end);   

columns = length(distributions_names);
x = [2, 3, 5];
linewidth = 1.5;
markersize = 1;

figure();
for it_columns=1:columns    
    pl = subplot(1, columns, it_columns);    
    bar(x, log10([data_america(:,it_columns), data_sa4ica(:,it_columns), data_GLICA(:,it_columns), data_QICA(:,it_columns)]),'hist');        
    if(it_columns==1)
        legend('america','sa4ica','GLICA','QICA','location','northwest');         
    end
    title(distributions_names{it_columns});
    grid on;
    xlabel('P');
    if(it_columns==1)
        ylabel(sprintf('Runtime [log_{10}(s)]\nK=%d, T=10^{%d}',n_sources(end), log10(n_samples(end))));
    end
end
