clear;
close all;
addpath('sim_data/');





load('sim_data_start_2022_02_17_10:51_end_2022_02_17_11:28_sim_total_corr_pure_ICA.mat');%P=2
data_america(1,:) = median(america_trial_time(:,1,3,end,:),5);   
data_QICA(1,:) = median(QICA_trial_time(:,1,3,end,:),5);
data_QICA_ex(1,:) = median(QICA_ex_trial_time(:,1,3,end,:),5);
data_sa4ica(1,:) = median(sa4ica_trial_time(:,1,3,end,:),5);   

load('sim_data_start_2022_02_17_12:01_end_2022_02_17_20:49_sim_total_corr_pure_ICA.mat');%P=3
data_america(2,:) = median(america_trial_time(:,1,3,end,:),5);   
data_QICA(2,:) = median(QICA_trial_time(:,1,3,end,:),5);
data_QICA_ex(2,:) = median(QICA_ex_trial_time(:,1,3,end,:),5);
data_sa4ica(2,:) = median(sa4ica_trial_time(:,1,3,end,:),5);    

load('sim_data_start_2022_02_18_10:14_end_2022_02_18_23:44_sim_total_corr_pure_ICA.mat');%P=5
data_america(3,:) = median(america_trial_time(:,1,3,end,:),5);   
data_QICA(3,:) = median(QICA_trial_time(:,1,3,end,:),5);
data_QICA_ex(3,:) = median(QICA_ex_trial_time(:,1,3,end,:),5);
data_sa4ica(3,:) = median(sa4ica_trial_time(:,1,3,end,:),5);     

columns = length(distributions_names);
x = [2, 3, 5];
linewidth = 1.5;

figure();
for it_columns=1:columns    
    pl = subplot(1, columns, it_columns);    
    semilogy(x, data_america(:,it_columns), '-*', ...
                x, data_sa4ica(:,it_columns), '--o', ...
                x, data_QICA(:,it_columns), '-.s', ...                
                x, data_QICA_ex(:,it_columns), '-.^', ...                
                'LineWidth', linewidth);
    if(it_columns==1)
        legend('america/GLICA','sa4ica','QICA','QICA-Exhaustive','location','northwest');         
    end
    title(distributions_names{it_columns});
    grid on;
    xlabel('P');
    if(it_columns==1)
        ylabel(sprintf('Median Runtime [s]\nK=%d, T=10^{%d}',n_sources(3), log10(n_samples(end))));
    end
end
