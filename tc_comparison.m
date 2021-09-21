clear;
close all;

ranking_america = [];
ranking_GLICA = [];
ranking_QICA = [];
ranking_sa4ica = [];
ranking_order = [];
ranking_QICA_ex = [];

load('sim_data_start_2021_08_31_09:29_end_2021_08_31_09:49_sim_total_corr_pure_ICA.mat');%P=2
lines = length(n_sources);
columns = length(distributions_names);
for it_lines=1:lines
    for it_columns=1:columns
        data_america = reshape(america_mean_total_corr_results(it_columns,:,it_lines,:),1,length(n_samples));   
        data_GLICA = reshape(GLICA_mean_total_corr_results(it_columns,:,it_lines,:),1,length(n_samples));   
        data_QICA = reshape(QICA_mean_total_corr_results(it_columns,:,it_lines,:),1,length(n_samples));   
        data_QICA_ex = reshape(QICA_ex_mean_total_corr_results(it_columns,:,it_lines,:),1,length(n_samples));   
        data_sa4ica = reshape(sa4ica_mean_total_corr_results(it_columns,:,it_lines,:),1,length(n_samples));   
        data_order = reshape(order_mean_total_corr_results(it_columns,:,it_lines,:),1,length(n_samples));
        
        B = [data_america; data_GLICA; data_QICA; data_sa4ica; data_order; data_QICA_ex];
        [V, R] = sort(B);
        Ties = B==min(B);
        keep_columns = sum(Ties)<=2; 
        R = R(:,keep_columns);
        Places = repmat((1:size(R,1))',1,size(R,2));
        
        ranking_america = [ranking_america; Places(R==1)];
        ranking_GLICA = [ranking_GLICA; Places(R==2)];
        ranking_QICA = [ranking_QICA; Places(R==3)];
        ranking_sa4ica = [ranking_sa4ica; Places(R==4)];
        ranking_order = [ranking_order; Places(R==5)];
        ranking_QICA_ex = [ranking_QICA_ex; Places(R==6)];
    end
end

fprintf('P=2 - Results - %d experiments:\n', length(ranking_america));
fprintf('AMERICA:\t%.2f (%.2f)\n', mean(ranking_america), std(ranking_america));
fprintf('GLICA:\t\t%.2f (%.2f)\n', mean(ranking_GLICA), std(ranking_GLICA));
fprintf('QICA:\t\t%.2f (%.2f)\n', mean(ranking_QICA), std(ranking_QICA));
fprintf('QICA-Exhaus.:\t%.2f (%.2f)\n', mean(ranking_QICA_ex), std(ranking_QICA_ex));
fprintf('SA4ICA:\t\t%.2f (%.2f)\n', mean(ranking_sa4ica), std(ranking_sa4ica));
fprintf('Order:\t\t%.2f (%.2f)\n', mean(ranking_order), std(ranking_order));

ranking_america = [];
ranking_GLICA = [];
ranking_QICA = [];
ranking_sa4ica = [];
ranking_QICA_ex = [];
load('sim_data_start_2021_08_31_16:37_end_2021_08_31_20:27_sim_total_corr_pure_ICA.mat');%P=3
lines = length(n_sources);
columns = length(distributions_names);
for it_lines=1:lines
    for it_columns=1:columns
        data_america = reshape(america_mean_total_corr_results(it_columns,:,it_lines,:),1,length(n_samples));   
        data_GLICA = reshape(GLICA_mean_total_corr_results(it_columns,:,it_lines,:),1,length(n_samples));   
        data_QICA = reshape(QICA_mean_total_corr_results(it_columns,:,it_lines,:),1,length(n_samples));   
        data_QICA_ex = reshape(QICA_ex_mean_total_corr_results(it_columns,:,it_lines,:),1,length(n_samples));   
        data_sa4ica = reshape(sa4ica_mean_total_corr_results(it_columns,:,it_lines,:),1,length(n_samples));           
        
        B = [data_america; data_GLICA; data_QICA; data_sa4ica; data_QICA_ex];
        [V, R] = sort(B);
        Ties = B==min(B);
        keep_columns = sum(Ties)<=2;
        R = R(:,keep_columns);
        Places = repmat((1:size(R,1))',1,size(R,2));
        
        ranking_america = [ranking_america; Places(R==1)];
        ranking_GLICA = [ranking_GLICA; Places(R==2)];
        ranking_QICA = [ranking_QICA; Places(R==3)];
        ranking_sa4ica = [ranking_sa4ica; Places(R==4)];
        ranking_QICA_ex = [ranking_QICA_ex; Places(R==5)];
    end
end

fprintf('P=3 - Results - %d experiments:\n', length(ranking_america));
fprintf('AMERICA:\t%.2f (%.2f)\n', mean(ranking_america), std(ranking_america));
fprintf('GLICA:\t\t%.2f (%.2f)\n', mean(ranking_GLICA), std(ranking_GLICA));
fprintf('QICA:\t\t%.2f (%.2f)\n', mean(ranking_QICA), std(ranking_QICA));
fprintf('QICA-Exhaus.:\t%.2f (%.2f)\n', mean(ranking_QICA_ex), std(ranking_QICA_ex));
fprintf('SA4ICA:\t\t%.2f (%.2f)\n', mean(ranking_sa4ica), std(ranking_sa4ica));

load('sim_data_start_2021_09_01_09:50_end_2021_09_01_18:41_sim_total_corr_pure_ICA.mat');%P=5
ranking_america = [];
ranking_GLICA = [];
ranking_QICA = [];
ranking_sa4ica = [];
lines = length(n_sources);
columns = length(distributions_names);
for it_lines=1:lines
    for it_columns=1:columns
        data_america = reshape(america_mean_total_corr_results(it_columns,:,it_lines,:),1,length(n_samples));   
        data_GLICA = reshape(GLICA_mean_total_corr_results(it_columns,:,it_lines,:),1,length(n_samples));   
        data_QICA = reshape(QICA_mean_total_corr_results(it_columns,:,it_lines,:),1,length(n_samples));   
        data_sa4ica = reshape(sa4ica_mean_total_corr_results(it_columns,:,it_lines,:),1,length(n_samples));                   
        
        B = [data_america; data_GLICA; data_QICA; data_sa4ica];
        [V, R] = sort(B);
        Ties = B==min(B);
        keep_columns = sum(Ties)<=2;
        R = R(:,keep_columns);
        Places = repmat((1:size(R,1))',1,size(R,2));
        
        ranking_america = [ranking_america; Places(R==1)];
        ranking_GLICA = [ranking_GLICA; Places(R==2)];
        ranking_QICA = [ranking_QICA; Places(R==3)];
        ranking_sa4ica = [ranking_sa4ica; Places(R==4)];
    end
end

fprintf('P=5 - Results - %d experiments:\n', length(ranking_america));
fprintf('AMERICA:\t%.2f (%.2f)\n', mean(ranking_america), std(ranking_america));
fprintf('GLICA:\t\t%.2f (%.2f)\n', mean(ranking_GLICA), std(ranking_GLICA));
fprintf('QICA:\t\t%.2f (%.2f)\n', mean(ranking_QICA), std(ranking_QICA));
fprintf('SA4ICA:\t\t%.2f (%.2f)\n', mean(ranking_sa4ica), std(ranking_sa4ica));