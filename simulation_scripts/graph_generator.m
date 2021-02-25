data = reshape(mean_total_corr_results,length(algorithms_names),length(n_samples));
x = log2(n_samples);
plot(x,data(1,:),x,data(2,:),x,data(3,:),x,data(4,:),x,data(5,:));legend('america','sa4ica','QICA','GLICA','order');
