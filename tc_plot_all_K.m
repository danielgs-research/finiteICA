clear;
close all;
addpath('sim_data/');


% load('sim_data_start_2021_08_31_16:37_end_2021_08_31_20:27_sim_total_corr_pure_ICA.mat');%P=3
% load('sim_data_start_2021_09_01_09:50_end_2021_09_01_18:41_sim_total_corr_pure_ICA.mat');%P=5
load('sim_data_start_2021_08_31_09:29_end_2021_08_31_09:49_sim_total_corr_pure_ICA.mat');%P=2


% components = [1 floor((1+length(n_sources))/2) length(n_sources)];
lines = length(n_sources);
columns = length(distributions_names);
x = n_samples;
linewidth = 2;
P = some_primes(1);

figure();
for it_lines=1:lines
    K = n_sources(it_lines);
    for it_columns=1:columns
        data_america = reshape(america_mean_total_corr_results(it_columns,:,it_lines,:),1,length(n_samples));   

        data_GLICA = reshape(GLICA_mean_total_corr_results(it_columns,:,it_lines,:),1,length(n_samples));   

        data_QICA = reshape(QICA_mean_total_corr_results(it_columns,:,it_lines,:),1,length(n_samples));   

        data_QICA_ex = reshape(QICA_ex_mean_total_corr_results(it_columns,:,it_lines,:),1,length(n_samples));   

        data_sa4ica = reshape(sa4ica_mean_total_corr_results(it_columns,:,it_lines,:),1,length(n_samples));   

        data_order = reshape(order_mean_total_corr_results(it_columns,:,it_lines,:),1,length(n_samples));

        pl = subplot(lines, columns,(it_lines-1)*columns + it_columns);
        if(P==2)
            semilogx(x, data_america, '-*', ...
                    x, data_sa4ica, '--o', ...
                    x, data_GLICA, ':d', ... 
                    x, data_QICA, '-.s', ...
                    x, data_QICA_ex, '--^', ...
                    x, data_order, '-x', 'LineWidth', linewidth);                
            if(it_lines==1 && it_columns==1)
                legend('america','sa4ica','GLICA','QICA','QICA-Exhaustive','order');            
            end
        else
            if(K<=4 && P==3)
                semilogx(x, data_america, '-*', ...
                        x, data_sa4ica, '--o', ...
                        x, data_GLICA, ':d', ... 
                        x, data_QICA_ex, '--^', ...
                        x, data_QICA, '-.s', 'LineWidth', linewidth);                
                if(it_columns==1)
                    legend('america','sa4ica','GLICA','QICA','QICA-Exhaustive');           
                end
            else
                semilogx(x, data_america, '-*', ...
                        x, data_sa4ica, '--o', ...
                        x, data_GLICA, ':d', ...                         
                        x, data_QICA, '-.s', 'LineWidth', linewidth);                
                if(it_columns==1)
                    legend('america','sa4ica','GLICA','QICA');           
                end
            end
        end
        if(it_lines==1)
            title(distributions_names{it_columns});            
        end            
        if(it_lines==lines)
            xlabel('T');
        end
        if(it_columns==1)
            ylabel(sprintf('Total Correlation [bits]\nK=%d', K));
        end
        grid on;
    end
end