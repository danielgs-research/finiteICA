function [stddata, meandata, data] = summary_statistics(data, space, round_digits)

data = reshape(data,space);
meandata = round(mean(data, length(space)),round_digits);
stddata = round(std(data, 0, length(space)),round_digits);


end