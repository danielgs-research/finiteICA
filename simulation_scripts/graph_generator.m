clear;
close all;
addpath('sim_data/');
load('sim_data_start_2021_07_09_15:48_end_2021_07_09_15:54_sim_total_corr_pure_ICA.mat');

lines = length(n_sources);
columns = length(distributions_names);
x = n_samples;
linewidth = 1.5;
markersize = 1;

figure();
for it_lines=1:lines
    for it_columns=1:columns
        data_america = reshape(america_mean_total_corr_results(it_columns,:,it_lines,:),1,length(n_samples));   
        data_GLICA = reshape(GLICA_mean_total_corr_results(it_columns,:,it_lines,:),1,length(n_samples));   
        data_QICA = reshape(QICA_mean_total_corr_results(it_columns,:,it_lines,:),1,length(n_samples));   
        data_sa4ica = reshape(sa4ica_mean_total_corr_results(it_columns,:,it_lines,:),1,length(n_samples));   
        data_order = reshape(order_mean_total_corr_results(it_columns,:,it_lines,:),1,length(n_samples));
        pl = subplot(lines, columns,(it_lines-1)*columns + it_columns);
        semilogx(x, data_america, '-*', 'LineWidth', linewidth, ...
                x, data_sa4ica, '--o', 'LineWidth', linewidth, ...
                x, data_GLICA, ':d', 'LineWidth', linewidth, ... 
                x, data_QICA, '-.s', 'LineWidth', linewidth, ...
                x, data_order, '-x', 'LineWidth', linewidth);
        legend('america','sa4ica','GLICA','QICA','order');        
        if(it_lines==1)
            title(sprintf('%s\nK=%d', distributions_names{it_columns},n_sources(it_lines)));
        else
            title(sprintf('K=%d', n_sources(it_lines)));
        end
%         xlabel('log_{10}(T)');
        xlabel('T');
        ylabel('Total Correlation [bits]');
    end
end