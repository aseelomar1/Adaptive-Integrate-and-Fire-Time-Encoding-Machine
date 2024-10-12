% genSinc generates a signal that is the sum of 2N+1 sinc functions.
% 
% Syntax:
%   [sig, td, cn] = genSinc(f, tds, TP, N)
%
% Inputs:
%   - f : Central frequency in Hz (scalar)
%   - tds : Time resolution (sampling interval) for the time vector (scalar)
%   - TP : Total time period for the time vector (scalar)
%   - N : Number of sinc functions on each side of the center (scalar)
%
% Outputs:
%   - sig : Generated signal (vector)
%   - td : Time vector (from 0 to TP with step tds) (vector)
%   - cn : Random amplitudes for the sinc components (vector)
%
% Description:
% This function generates a signal that is a sum of 2N+1 sinc functions. 
% The signal is constructed as:
%   sig = sum(cn(n) * sinc((td - n * ts) / ts))
% where ts is the sampling interval determined by the frequency f, and cn 
% are random amplitudes. The output includes the generated signal, the time 
% vector, and the vector of random amplitudes.
%
% Example usage:
%   f = 30;
%   tds = 0.001;
%   TP = 1;
%   N = 10;
%   [sig, td, cn] = genSinc(f, tds, TP, N);
%
% Author: [Aseel Omar]

function [sig, td, cn] = genSinc(f, tds, TP, N)

    % Sampling frequency and interval
    Fs = 2 * f;           % Sampling frequency (Nyquist rate)
    ts = 1 / Fs;          % Sampling interval

    % Time vector from 0 to TP with step tds
    td = (0:tds:TP).';    % Time vector (column vector)

    % Initialize signal vector to zeros
    sig = zeros(length(td), 1);

    % Generate random amplitudes for the sinc components
    cn = rand(2 * N + 1, 1);  % Random amplitudes

    % Sum 2N+1 sinc functions to generate the signal
    for n = -N:N
        % Add each sinc component multiplied by cn(n + N + 1) to the signal
        sig = sig + cn(n + N + 1) * sinc((td - n * ts) / ts);
    end
end
