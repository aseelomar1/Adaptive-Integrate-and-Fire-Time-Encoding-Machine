% compute_G_matrix Compute generating matrix for TEM simulation
%
% Syntax:
%   [G, thetan] = compute_G_matrix(t_vec, tn, bw)
%
% Inputs:
%   - t_vec: Time grid (vector)
%   - tn: Spike times (vector)
%   - bw: Maximal frequency (bandwidth) in radians/second (scalar)
%
% Outputs:
%   - G: Generating matrix for signal recovery (matrix)
%   - thetan: Midpoints between consecutive spike times (vector)
%
% Description:
%   This function computes the generating matrix `G` for the time encoding
%   machine (TEM) simulation, along with the midpoints `thetan` between consecutive
%   spike times `tn`. The generating matrix is constructed using the sinc
%   function based on the bandwidth 'bw' and the time grid `t_vec`.
%
% Example usage:
%   [G, thetan] = compute_G_matrix(t_vec, tn, OmegaMax);
%
% Author: [Aseel Omar]

function [G, thetan] = compute_G_matrix(t_vec, tn, bw)

    % Number of intervals between spike times
    sampleLen = length(tn) - 1;
    
    % Initialize thetan as the midpoints between consecutive spike times
    thetan = zeros(sampleLen, 1);
    for kk = 1:sampleLen
        thetan(kk) = (tn(kk) + tn(kk+1)) / 2;
    end
    
    % Number of samples in thetan
    samples = numel(thetan);
    
    % Time step size (assuming uniform spacing in t_vec)
    dt = t_vec(2) - t_vec(1);
    
    % Define the scaled sinc function
    MySinc = @(x) (bw/pi) * sinc((bw * x) / pi);
    
    % Initialize generating matrix G
    G = zeros(sampleLen, samples);
    
    % Compute the generating matrix G entries
    for nn = 1:sampleLen
        % Select the time points in t_vec between tn(nn) and tn(nn+1)
        t_used = t_vec((t_vec > tn(nn)) & (t_vec <= tn(nn+1)));
        
        % Ensure t_used is a column vector
        if size(t_used, 1) == 1
            t_used = t_used';
        end
        
        % Create a matrix where each row is t_used minus thetan
        t_used_mat = repmat(t_used, 1, samples)';
        shiftSinc = MySinc(t_used_mat - thetan);  % Compute the sinc values
        
        % Sum over time points and multiply by dt to get the matrix entry
        G(nn, :) = dt * sum(shiftSinc, 2);
    end

end
