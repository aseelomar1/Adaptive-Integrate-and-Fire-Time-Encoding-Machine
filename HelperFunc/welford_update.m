function [new_mean, new_var] = welford_update(prev_mean, prev_var, n, new_data)
    % Welford's method for online variance calculation (population variance).
    % This function updates the mean and variance given new data.
    %
    % Input:
    %   prev_mean: Mean of data up to the (n-1)th data point
    %   prev_var: Variance of data up to the (n-1)th data point
    %   n: Number of data points seen so far (including new_data)
    %   new_data: New data point to incorporate
    %
    % Output:
    %   new_mean: Updated mean
    %   new_var: Updated variance

    delta = new_data - prev_mean;   % Difference between new data and current mean
    new_mean = prev_mean + delta / n;   % Update mean
    delta2 = new_data - new_mean;   % Recalculate the difference using new mean
    new_var = ((n - 1) * prev_var + delta * delta2) / n;   % Update variance
end
