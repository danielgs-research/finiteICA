clear;
close all;
addpath('sim_data/');





load('sim_data_start_2021_08_31_09:29_end_2021_08_31_09:49_sim_total_corr_pure_ICA.mat');%P=2
data_america(1,:) = median(america_trial_time(:,1,3,end,:),5);   
data_GLICA(1,:) = GLICA_mean_trial_time(:,1,3,end);   
data_QICA(1,:) = QICA_mean_trial_time(:,1,3,end);
data_QICA_ex(1,:) = QICA_ex_mean_trial_time(:,1,3,end);
data_sa4ica(1,:) = sa4ica_mean_trial_time(:,1,3,end);   

load('sim_data_start_2021_08_31_16:37_end_2021_08_31_20:27_sim_total_corr_pure_ICA.mat');%P=3
data_america(2,:) = median(america_trial_time(:,1,3,end,:),5);   
data_GLICA(2,:) = GLICA_mean_trial_time(:,1,3,end);   
data_QICA(2,:) = QICA_mean_trial_time(:,1,3,end);   
data_QICA_ex(2,:) = QICA_ex_mean_trial_time(:,1,3,end);
data_sa4ica(2,:) = sa4ica_mean_trial_time(:,1,3,end);   

load('sim_data_start_2021_09_01_09:50_end_2021_09_01_18:41_sim_total_corr_pure_ICA.mat');%P=5
data_america(3,:) = america_mean_trial_time(:,1,3,end);   
data_GLICA(3,:) = GLICA_mean_trial_time(:,1,3,end);   
data_QICA(3,:) = QICA_mean_trial_time(:,1,3,end);   
data_QICA_ex(3,:) = QICA_ex_mean_trial_time(:,1,3,end);
data_sa4ica(3,:) = sa4ica_mean_trial_time(:,1,3,end);   

columns = length(distributions_names);
x = [2, 3, 5];
linewidth = 1.5;

figure();
for it_columns=1:columns    
    pl = subplot(1, columns, it_columns);    
    semilogy(x, data_america(:,it_columns), '-*', ...
                x, data_sa4ica(:,it_columns), '--o', ...
                x, data_GLICA(:,it_columns), ':d', ... 
                x, data_QICA(:,it_columns), '-.s', ...                
                x, data_QICA_ex(:,it_columns), '-.^', ...                
                'LineWidth', linewidth);
    if(it_columns==1)
        legend('america','sa4ica','GLICA','QICA','QICA-Exhaustive','location','northwest');         
    end
    title(distributions_names{it_columns});
    grid on;
    xlabel('P');
    if(it_columns==1)
        ylabel(sprintf('Median Runtime [s]\nK=%d, T=10^{%d}',n_sources(3), log10(n_samples(end))));
    end
end
