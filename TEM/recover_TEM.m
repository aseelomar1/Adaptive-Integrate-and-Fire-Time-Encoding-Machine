% recover_TEM Decode a signal encoded with an IF-TEM or AIF-TEM.
%
% Syntax:
%   recSig = recover_TEM(td, tn, bw, b, kappa, th)
%
% Inputs:
%   - td: Time vector (vector)
%   - tn: Trigger times from A/IF-TEM output (vector)
%   - bw: Bandwidth of the signal in [rad/s] (scalar)
%   - b: Bias for the A/IF-TEM (vector or scalar)
%   - kappa: Integrator acceleration parameter (scalar)
%   - th: Threshold (scalar)
%
% Outputs:
%   - recSig: Recovered signal (vector)
%
% Description:
%   This function decodes a signal encoded with an Adaptive or Integrate-and-Fire
%   Time Encoding Machine (A/IF-TEM). It uses the time vector `td` and trigger
%   times `tn`, along with bandwidth `bw`, bias `b`, and the A/IF-TEM parameters
%   `kappa` and `th` (threshold). The recovered signal is returned as `recSig`.
%
% Example usage:
%   recSig = recover_TEM(td, tn, bw, b, kappa, th);
%
% Author: [Aseel Omar]

function recSig = recover_TEM(td, tn, bw, b, kappa, th)

    % Calculate the difference between consecutive trigger times (intervals)
    tn_intervals = tn(2:end) - tn(1:end-1);
    
    % Compute the generating matrix and thetan (midpoints between spikes)
    [G, thetan] = compute_G_matrix(td, tn, bw);
    
    % Ensure the bias vector b is the correct orientation
    if size(tn_intervals, 1) == size(b, 1) && size(tn_intervals, 2) == size(b, 2) && numel(b) > 1
        b = b'; % Transpose bias if needed
    end
    
    % Calculate the reconstruction coefficients (q)
    q = (kappa * th - b' .* tn_intervals);
    
    % Initialize the recovered signal as a zero vector
    recSig = zeros(size(td));
    
    % Compute the pseudoinverse of the recovery matrix
    G_inv = pinv(G);
    
    %
    if(size(q,1) ~= size(G_inv,1))
        q=q';
    end

    % Calculate the coefficient vector
    c_coeff = G_inv * q;
    
    % Number of sinc terms for reconstruction
    nsh = length(thetan);
    
    % Bandwidth divided by pi
    bwpi = bw / pi;
    
    % Reconstruct the signal by summing weighted sinc functions
    for i = 1:nsh
        recSig = recSig + sinc(bwpi * (td - thetan(i))) * bwpi * c_coeff(i);
    end

end
