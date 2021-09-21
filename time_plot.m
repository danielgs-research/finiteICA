clear;
close all;
addpath('sim_data/');

load('sim_data_start_2021_09_01_09:50_end_2021_09_01_18:41_sim_total_corr_pure_ICA.mat');%P=5
% load('sim_data_start_2021_08_31_09:29_end_2021_08_31_09:49_sim_total_corr_pure_ICA.mat');%P=2
% load('sim_data_start_2021_08_31_16:37_end_2021_08_31_20:27_sim_total_corr_pure_ICA.mat');%P=3



columns = length(distributions_names);
x = n_sources;
linewidth = 1.5;
P = some_primes(1);

figure();

for it_columns=1:columns
    data_america = reshape(median(america_trial_time(it_columns,:,:,end,:),5),1,length(n_sources));   
    data_GLICA = reshape(median(GLICA_trial_time(it_columns,:,:,end,:),5),1,length(n_sources));   
    data_QICA = reshape(median(QICA_trial_time(it_columns,:,:,end,:),5),1,length(n_sources));   
    data_sa4ica = reshape(median(sa4ica_mean_trial_time(it_columns,:,:,end,:),5),1,length(n_sources));       
    data_order = reshape(median(order_trial_time(it_columns,:,:,end,:),5),1,length(n_sources));
    data_QICA_ex = reshape(median(QICA_ex_trial_time(it_columns,:,:,end,:),5),1,length(n_sources));   
    pl = subplot(1, columns, it_columns);
    if(P==2)    
        semilogy(x, data_america, '-*', ...
                    x, data_sa4ica, '--o', ...
                    x, data_GLICA, ':d', ... 
                    x, data_QICA, '-.s', ...
                    x, data_QICA_ex, '--^', ...
                    x, data_order, '-x', 'LineWidth', linewidth);
        if(it_columns==1)
            legend('america','sa4ica','GLICA','QICA','QICA-Exhaustive','order','location','northwest');         
        end
    else
        if(P==3)
            semilogy(x, data_america, '-*', ...
                    x, data_sa4ica, '--o', ...
                    x, data_GLICA, ':d', ... 
                    x, data_QICA, '-.s', ...
                    x, data_QICA_ex, '--^', ...
                    'LineWidth', linewidth);
            if(it_columns==1)
                legend('america','sa4ica','GLICA','QICA','QICA-Exhaustive','location','northwest');         
            end
        else
            semilogy(x, data_america, '-*', ...
                        x, data_sa4ica, '--o', ...
                        x, data_GLICA, ':d', ... 
                        x, data_QICA, '-.s', ...
                        'LineWidth', linewidth);        
            if(it_columns==1)
                legend('america','sa4ica','GLICA','QICA','location','northwest');         
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
