function b_q = bias_quant(b, Delta,bmin,bmax)
    % Quantizes the input bias b based on the step size Delta.
    % Input:
    %   b: Bias value to be quantized
    %   Delta: Quantization step size
    % Output:
    %   b_q: Quantized bias value
    
    k = round(b / Delta + 0.5);  % Round to nearest quantized step
    b_q = k * Delta+bmin;  % Quantized bias
    b_q = min(b_q,bmax);
end
