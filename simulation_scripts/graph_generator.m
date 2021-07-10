clear;
close all;
addpath('sim_data/');
load('sim_data_start_2021_07_09_15:48_end_2021_07_09_15:54_sim_total_corr_pure_ICA.mat');

lines = length(n_sources);
columns = length(distributions_names);
% x = log10(n_samples);
x = n_samples;

figure();
for it_lines=1:lines
    for it_columns=1:columns
        data_america = reshape(america_mean_total_corr_results(it_columns,:,it_lines,:),1,length(n_samples));   
        data_GLICA = reshape(GLICA_mean_total_corr_results(it_columns,:,it_lines,:),1,length(n_samples));   
        data_QICA = reshape(QICA_mean_total_corr_results(it_columns,:,it_lines,:),1,length(n_samples));   
        data_sa4ica = reshape(sa4ica_mean_total_corr_results(it_columns,:,it_lines,:),1,length(n_samples));   
        data_order = reshape(order_mean_total_corr_results(it_columns,:,it_lines,:),1,length(n_samples));
        pl = subplot(lines, columns,(it_lines-1)*columns + it_columns);
        semilogx(x, data_america, '-*', x, data_sa4ica, '--o', ...
                x, data_GLICA, ':d', x, data_QICA, '-.s', x, data_order, '-x');
        legend('america','sa4ica','GLICA','QICA','order');
        title(sprintf('%s - K=%d', string(distributions_names(it_columns)),n_sources(it_lines)));
%         xlabel('log_{10}(T)');
        xlabel('T');
        ylabel('Total Correlation [bits]');
    end
end