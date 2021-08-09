clear;
close all;
addpath('sim_data/');
load('sim_data_start_2021_07_10_12:29_end_2021_07_12_15:57_sim_total_corr_pure_ICA.mat');%P=5
% load('sim_data_start_2021_07_09_22:13_end_2021_07_10_12:25_sim_total_corr_pure_ICA.mat');%P=3
% load('sim_data_start_2021_07_09_15:48_end_2021_07_09_15:54_sim_total_corr_pure_ICA.mat');%P=2

lines = length(n_sources);
columns = length(distributions_names);
winners = [];
for it_lines=1:lines
    for it_columns=1:columns
        data_america = reshape(america_mean_total_corr_results(it_columns,:,it_lines,:),1,length(n_samples));   
        data_GLICA = reshape(GLICA_mean_total_corr_results(it_columns,:,it_lines,:),1,length(n_samples));   
        data_QICA = reshape(QICA_mean_total_corr_results(it_columns,:,it_lines,:),1,length(n_samples));   
        data_sa4ica = reshape(sa4ica_mean_total_corr_results(it_columns,:,it_lines,:),1,length(n_samples));   
        data_order = reshape(order_mean_total_corr_results(it_columns,:,it_lines,:),1,length(n_samples));
        
        B = [data_america; data_GLICA; data_QICA; data_sa4ica; data_order];
        [V, R] = sort(B);
        
        disp(V);
        disp(R);
        
%         winners = [winners w(sum(T)==1)];        
%         disp(B)
%         disp(T)
%         disp(sum(T))
%         disp(winners)
        pause
    end
end

fprintf('AMERICA: %d wins\n', sum(winners==1));
% fprintf('GLICA: %d wins\n', sum(winners==2));
% fprintf('QICA: %d wins\n', sum(winners==3));
% fprintf('SA4ICA: %d wins\n', sum(winners==4));
% fprintf('Order: %d wins\n', sum(winners==5));

