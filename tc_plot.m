clear;
close all;
addpath('sim_data/');

% load('sim_data_start_2021_08_06_08:19_end_2021_08_06_08:29_sim_total_corr_pure_ICA.mat');%P=2
% load('sim_data_start_2021_08_06_08:40_end_2021_08_06_09:03_sim_total_corr_pure_ICA.mat');%P=3
load('sim_data_start_2021_08_06_12:41_end_2021_08_06_21:27_sim_total_corr_pure_ICA.mat');%P=5

lines = length(n_sources);
columns = length(distributions_names);
x = n_samples;
linewidth = 2;
markersize = 1;
binaryP = some_primes(1)==2;

figure();
for it_lines=1:lines
    for it_columns=1:columns
        data_america = reshape(america_mean_total_corr_results(it_columns,:,it_lines,:),1,length(n_samples));   
        data_GLICA = reshape(GLICA_mean_total_corr_results(it_columns,:,it_lines,:),1,length(n_samples));   
        data_QICA = reshape(QICA_mean_total_corr_results(it_columns,:,it_lines,:),1,length(n_samples));   
        data_sa4ica = reshape(sa4ica_mean_total_corr_results(it_columns,:,it_lines,:),1,length(n_samples));   
        data_order = reshape(order_mean_total_corr_results(it_columns,:,it_lines,:),1,length(n_samples));
        pl = subplot(lines, columns,(it_lines-1)*columns + it_columns);
        if(binaryP)
            semilogx(x, data_america, '-*', ...
                    x, data_sa4ica, '--o', ...
                    x, data_GLICA, ':d', ... 
                    x, data_QICA, '-.s', ...
                    x, data_order, '-x', 'LineWidth', linewidth);                
            if(it_lines==1 && it_columns==1)
                legend('america','sa4ica','GLICA','QICA','order');            
            end
        else
            semilogx(x, data_america, '-*', ...
                    x, data_sa4ica, '--o', ...
                    x, data_GLICA, ':d', ... 
                    x, data_QICA, '-.s', 'LineWidth', linewidth);                
            if(it_lines==1 && it_columns==1)
                legend('america','sa4ica','GLICA','QICA');           
            end
        end
        if(it_lines==1)
            title(distributions_names{it_columns});            
        end            
        if(it_lines==lines)
            xlabel('T');
        end
        if(it_columns==1)
            ylabel(sprintf('Total Correlation [bits]\nK=%d',n_sources(it_lines)));
        end
        grid on;
    end
end