function h = entropy_from_frequencies(ffq, log_base)
	% FFQ must be already in frequency form, integer values are not counted


	% IMPORTANT: THE FREQUENCIES MUST BE COLUMN-WISE
	% I.E assert(all(sum(ffq,1) == 1))


	% since it's double values and we have float problems,
	% we will assert the way below.
    
    if(nargin<2) %MATLAB
        log_base = 2;
    end

	assert(all(sum(ffq,1) >= 0.95));
	assert(all(sum(ffq,1) <= 1.1));


	h=-sum(ffq.*log2(ffq+eps),1)./log2(log_base);

end
