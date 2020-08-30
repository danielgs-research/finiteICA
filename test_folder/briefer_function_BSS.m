% First, load your saved data manually
% load('file_name_in_your_actual_path')

% I believe we can do this differently, without many levels
if (exist('some_primes'))
	disp(some_primes)
	if(exist('n_sources'))
		disp(n_sources)
		if(exist('n_samples'))
			disp(n_samples)
			if(exist('algorithms_names'))
				disp(algorithms_names)
				if(exist('n_trials'))
					disp(n_trials)
                end
            end
        end
    end
end


a_str = sprintf('Total runs of algorithms is %d',n_cases)
disp(a_str)
a_str = sprintf('Kullback-leibler divergence threshold is %f',threshold)
disp(a_str)
a_str = sprintf('Simulation full name is %s',saved_sim)
disp(a_str)