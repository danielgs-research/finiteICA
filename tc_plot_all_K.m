clear;
close all;
addpath('sim_data/');


% load('sim_data_start_2022_02_18_10:14_end_2022_02_18_23:44_sim_total_corr_pure_ICA.mat');%P=5
% load('sim_data_start_2022_02_17_12:01_end_2022_02_17_20:49_sim_total_corr_pure_ICA.mat');%P=3
load('sim_data_start_2022_02_17_10:51_end_2022_02_17_11:28_sim_total_corr_pure_ICA.mat'); %P=2


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

        data_QICA = reshape(QICA_mean_total_corr_results(it_columns,:,it_lines,:),1,length(n_samples));   

        data_QICA_ex = reshape(QICA_ex_mean_total_corr_results(it_columns,:,it_lines,:),1,length(n_samples));   

        data_sa4ica = reshape(sa4ica_mean_total_corr_results(it_columns,:,it_lines,:),1,length(n_samples));   

        data_order = reshape(order_mean_total_corr_results(it_columns,:,it_lines,:),1,length(n_samples));

        pl = subplot(lines, columns,(it_lines-1)*columns + it_columns);
        if(P==2)
            semilogx(x, data_america, '-*', ...
                    x, data_sa4ica, '--o', ...                    
                    x, data_QICA, '-.s', ...
                    x, data_QICA_ex, '--^', ...
                    x, data_order, '-x', 'LineWidth', linewidth);                
            if(it_lines==1 && it_columns==1)
                legend('america/GLICA','sa4ica','QICA','QICA-Exhaustive','order');            
            end
        else
            if(K<=4 && P==3)
                semilogx(x, data_america, '-*', ...
                        x, data_sa4ica, '--o', ...                        
                        x, data_QICA_ex, '--^', ...
                        x, data_QICA, '-.s', 'LineWidth', linewidth);                
                if(it_columns==1)
                    legend('america/GLICA','sa4ica','QICA','QICA-Exhaustive');           
                end
            else
                semilogx(x, data_america, '-*', ...
                        x, data_sa4ica, '--o', ...                        
                        x, data_QICA, '-.s', 'LineWidth', linewidth);                
                if(it_columns==1)
                    legend('america/GLICA','sa4ica','QICA');           
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