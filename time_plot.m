clear;
close all;
addpath('sim_data/');

% load('sim_data_start_2021_07_10_12:29_end_2021_07_12_15:57_sim_total_corr_pure_ICA.mat');%P=5
% load('sim_data_start_2021_07_09_22:13_end_2021_07_10_12:25_sim_total_corr_pure_ICA.mat');%P=3
% load('sim_data_start_2021_07_09_15:48_end_2021_07_09_15:54_sim_total_corr_pure_ICA.mat');%P=2
load('sim_data_start_2021_08_06_08:40_end_2021_08_06_09:03_sim_total_corr_pure_ICA.mat');



columns = length(distributions_names);
x = n_sources;
linewidth = 1.5;
markersize = 1;
binaryP = some_primes(1)==2;

figure();

for it_columns=1:columns
    data_america = reshape(america_mean_trial_time(it_columns,:,:,end),1,length(n_sources));   
    data_GLICA = reshape(GLICA_mean_trial_time(it_columns,:,:,end),1,length(n_sources));   
    data_QICA = reshape(QICA_mean_trial_time(it_columns,:,:,end),1,length(n_sources));   
    data_sa4ica = reshape(sa4ica_mean_trial_time(it_columns,:,:,end),1,length(n_sources));   
    data_order = reshape(order_mean_trial_time(it_columns,:,:,end),1,length(n_sources));
    pl = subplot(1, columns, it_columns);
    if(binaryP)
        bar(x, [data_america; data_sa4ica; data_GLICA; data_QICA; data_order]','hist');                                    
        if(it_columns==1)
            legend('america','sa4ica','GLICA','QICA','order','location','northwest');         
        end
    else
        bar(x, [data_america; data_sa4ica; data_GLICA; data_QICA]','hist');                                    
        if(it_columns==1)
            legend('america','sa4ica','GLICA','QICA','location','northwest');         
        end
    end
    title(distributions_names{it_columns});
    grid on;
    xlabel('K');
    if(it_columns==1)
        ylabel(sprintf('Runtime [s]\nT=10^{%d}',log10(n_samples(end))));
    end
end
