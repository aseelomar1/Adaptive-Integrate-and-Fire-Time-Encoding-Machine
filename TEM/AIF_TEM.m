%AIF_TEM encode a signal using Adaptive Integrate and Fire Time Encoding 
% Machine. 
% [tr, bias] = AIF_TEM(x, dt, b0, th, k, alpha1, alpha2, beta, Delta, bmin, w, debug_mode)
% encodes the signal x using Adaptive Integrate and Fire Time Encoding 
%
% Inputs:
%   - x: Signal (vector) to be encoded
%   - dt: Signal resolution (scalar)
%   - b0: Initial Bias (scalar)
%   - th: Threshold (scalar)
%   - k: Integrator accelerator (optional, default = 1)
%   - alpha1: Parameter for Max Amplitude Predictor (optional, default = 0.98)
%   - alpha2: Parameter for Max Amplitude Predictor (optional, default = 1)
%   - beta: Bias adaptation factor (optional, default = 0)
%   - Delta: Quantization step for bias (optional, default = 1e-3)
%   - bmin: Minimum bias allowed (optional, default = 0)
%   - w: Window size for maximum bias calculation (optional, default = 1)
%   - debug_mode: Debugging mode (optional, default = 0)
%
% Outputs:
%   - tr: Trigger time indices (vector)
%   - bias: Bias values at the trigger times (vector)
%   - (Optional) y_out: Integrator output at each time step (vector, if debug_mode is on)
%   - (Optional) c: Real local maximum amplitude (vector, if debug_mode is on)
%
% Example:
%   [tr, bias] = AIF_TEM(x, dt, b0, th, 1, 0.98, 1, 0, 1e-3, 0, 1, 0);
%
% Author: Aseel Omar

function [varargout] = AIF_TEM(x, dt, b0, th, varargin)

    % Default parameters
    k = 1;
    alpha1 = 0.98;
    alpha2 = 1;
    beta = 0;
    Delta = 10^-3;
    bmin = 0;
    w = 1;
    intg_o = [];  % Initialize empty array
    bias = b0;   % Initial bias value
    debug_mode = 0;
    integOut = 0;
    c = [];
    
    % Parse optional arguments
    if nargin >= 5, k = varargin{1}; end
    if nargin >= 6, alpha1 = varargin{2}; end
    if nargin >= 7, alpha2 = varargin{3}; end
    if nargin >= 8, beta = varargin{4}; end
    if nargin >= 9, Delta = varargin{5}; end
    if nargin >= 10, bmin = varargin{6}; end
    if nargin >= 11, w = varargin{7}; end
    if nargin >= 12, debug_mode = varargin{8}; end
    if nargin >= 13, error('Too many input arguments.'); end
    
    if isempty(x)
        if debug_mode  
            varargout = {[], bias, intg_o, c};
        else 
            varargout = {[], bias};
        end
        return
    end
    
    % Initialization
    Nx = length(x);
    j = 1;
    c_est(1) = b0 - beta;   % Initial local amplitude estimate
    mean_val = b0 - beta;   % Mean value for local amplitude
    var_val = 0;            % Variance of local amplitude
    bj = b0;                % Current bias
    tr = [];                % Initialize empty array for trigger times
    b = [];                 % Initialize empty array for bias values in window

     % Define the function to compute integrator output at each step
    integOutcalc = @(integOut, i, bj) integOut + dt * (bj + x(i)) / k;

    % Main loop
    for i = 1:Nx
        % Compute integrator output
        integOut = integOutcalc(integOut, i, bj);
        intg_o(i) = integOut;
        
        % Threshold firing condition
        if integOut >= th
            tr(j) = i; % Record trigger time index
            
            % Calculate time difference since last trigger
            if j == 1
                T = (i - 1) * dt;
            else
                T = (i - tr(j - 1)) * dt;
            end
            
            % Max Amplitude Predictor
            z = abs(-bias(j) + (k * th / T));
            c_est = (1 - alpha1) * c_est + alpha1 * z;
            [mean_val, var_val] = welford_update(mean_val, var_val, j + 1, c_est); % Update mean and variance
            c_fut = abs(c_est) + alpha2 * sqrt(var_val); % Estimate future max amplitude
            
            % Bias Adaptor
            b(j+1) = max(c_fut + beta, bmin); % Update bias
            if j+1 <= w
                bias(j+1) = max(b(1:j+1));
            else
                bias(j+1) = max(b(j+1 - (w - 1):j+1));
            end
            bias(j+1) = bias_quant(bias(j+1), Delta,bmin,b0); % Quantize bias
            bj = bias(j+1); % Update current bias

            % Debugging mode: Record real local max amplitude
            if debug_mode
                if j == 1
                    c(1) = max(abs(x(1:i)));
                else
                    c(j) = max(abs(x(tr(j-1)+1:i)));
                end
            end
            
            % Reset integrator and increment trigger counter
            j = j + 1;
            integOut = integOut - th;
        end
    end
    
    % Output based on debug mode
    if debug_mode
        varargout = {tr, bias, intg_o, c};
    else
        varargout = {tr, bias};
    end

       % Handle case where no triggers were recorded
    if isempty(tr)
        warning('No triggers were recorded; signal did not reach the threshold.');
    end 

end
