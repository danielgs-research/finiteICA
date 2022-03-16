clear;
close all;
addpath('sim_data/');

load('sim_data_start_2022_02_18_10:14_end_2022_02_18_23:44_sim_total_corr_pure_ICA.mat');%P=5
% load('sim_data_start_2022_02_17_12:01_end_2022_02_17_20:49_sim_total_corr_pure_ICA.mat');%P=3
% load('sim_data_start_2022_02_17_10:51_end_2022_02_17_11:28_sim_total_corr_pure_ICA.mat'); %P=2



columns = length(distributions_names);
x = n_sources;
linewidth = 1.5;
P = some_primes(1);

figure();

for it_columns=1:columns
    data_america = reshape(median(america_trial_time(it_columns,:,:,end,:),5),1,length(n_sources));   
    data_QICA = reshape(median(QICA_trial_time(it_columns,:,:,end,:),5),1,length(n_sources));   
    data_sa4ica = reshape(median(sa4ica_mean_trial_time(it_columns,:,:,end,:),5),1,length(n_sources));       
    data_order = reshape(median(order_trial_time(it_columns,:,:,end,:),5),1,length(n_sources));
    data_QICA_ex = reshape(median(QICA_ex_trial_time(it_columns,:,:,end,:),5),1,length(n_sources));   
    pl = subplot(1, columns, it_columns);
    if(P==2)    
        semilogy(x, data_america, '-*', ...
                    x, data_sa4ica, '--o', ...
                    x, data_QICA, '-.s', ...
                    x, data_QICA_ex, '--^', ...
                    x, data_order, '-x', 'LineWidth', linewidth);

        if(it_columns==1)
            legend('america/GLICA','sa4ica','QICA','QICA-Exhaustive','order','location','northwest');              
        end
    else
        if(P==3)
            semilogy(x, data_america, '-*', ...
                    x, data_sa4ica, '--o', ...
                    x, data_QICA, '-.s', ...
                    x, data_QICA_ex, '--^', ...
                    'LineWidth', linewidth);               
            if(it_columns==1)
                legend('america/GLICA','sa4ica','QICA','QICA-Exhaustive','location','northwest');
            end
        else
            semilogy(x, data_america, '-*', ...
                        x, data_sa4ica, '--o', ...
                        x, data_QICA, '-.s', ...
                        'LineWidth', linewidth);         
            if(it_columns==1)
                legend('america/GLICA','sa4ica','QICA','location','northwest');
            end
        end
    end
    title(distributions_names{it_columns});
    grid on;
    xlabel('K');
    if(it_columns==1)
        ylabel(sprintf('Median Runtime [s]\nT=10^{%d}',log10(n_samples(end))));
    end
end
