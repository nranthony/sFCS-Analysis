function [corr_norm, corr_sdev, lags, A] = cross_corr_subsets(ps_times, num_subsets)
%CROSS_CORR_SUBSETS Repeats ACF of subsets to generate variance
%   splits picosecond times (ps_times) in num_subsets to be correlated to
%   generate average and standard deviation for ACF

ps_len = length(ps_times);
subset_wid = floor(ps_len/num_subsets);

ps_sub = ps_times(1:10);
[corr, lags] = cross_corr(ps_sub, ps_sub, 1000000, 10000000000000, 10, 0);

A = zeros(length(corr), num_subsets);


% Ques: do ps_times need to start from zero?

i=2;
for i = 1:num_subsets
    start = (i - 1) * subset_wid + 1;
    stop = start + subset_wid - 1;
    ps_sub = ps_times(start:stop);
    ps_sub = ps_sub - ps_sub(1);
    [corr, lags] = cross_corr(ps_sub, ps_sub, 1000000, 10000000000000, 10, 0);
    A(:,i) = corr;
end

lags = lags .* 1E-6; % converts to microseconds
corr_norm = mean(A,2);
corr_sdev = std(A,0,2);
corr_norm = corr_norm - 1; % convert from F to delta F

end
