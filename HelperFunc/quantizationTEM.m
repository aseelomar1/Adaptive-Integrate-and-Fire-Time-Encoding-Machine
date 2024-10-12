

function tn_quant = quantizationTEM(tn, tmax, tmin, bits,t0)
% quantizationTEM performs quantization on the time intervals between
% spike times (tn) using a specified number of quantization bits.
%
% Inputs:
%   - tn: Spike times (vector)
%   - tmax: Maximum value for quantization
%   - tmin: Minimum value for quantization
%   - bits: Number of bits for quantization
%   - t0: start time
%
% Output:
%   - tn_quant: Quantized time intervals after quantization (vector)
%
% Description:
%   This function calculates the time intervals between consecutive spike
%   times, quantizes them using a uniform quantization method, and returns
%   the cumulative sum of the quantized intervals.
%

    % Calculate the time intervals between consecutive spike times
    tn_intervals = tn(2:end) - tn(1:end-1);
    tn_intervals = [tn(1)-t0 , tn_intervals];
    % Number of quantization levels based on the number of bits
    quantizelevels = 2^bits;
    
    % Step size for each quantization level
    step_jump = (tmax - tmin) / (quantizelevels - 1);
    
    % Define the partition and codebook for quantization
    partition = tmin + step_jump/2 : step_jump : tmax - step_jump/2;
    codebook  = tmin : step_jump : tmax;
    
    % Perform quantization
    [index, quantized] = quantiz(tn_intervals, partition, codebook);
    
    % Return the cumulative sum of the quantized intervals
    tn_quant = cumsum(quantized') + t0;
   
end
