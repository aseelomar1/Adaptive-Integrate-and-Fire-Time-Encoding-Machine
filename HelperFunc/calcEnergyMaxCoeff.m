function [E, cmax] = calcEnergyMaxCoeff(f, N, unit)
% calcEnergyMaxCoeff calculates the energy and the maximum coefficient 
% for a signal built from 2N+1 sinc functions, supporting both Hz and kHz.
%
% Inputs:
%   - f: Frequency (scalar)
%   - N: Number of sinc functions on each side of the center (scalar)
%   - unit: String specifying the frequency unit ('Hz' or 'kHz')
%
% Outputs:
%   - E: Energy of the generated signal (scalar)
%   - cmax: Maximum coefficient related to the energy and frequency (scalar)
%
% Description:
%   The function generates a signal that is the sum of 2N+1 sinc functions.
%   It computes the energy of the signal using the sum of squared magnitudes
%   and returns the maximum coefficient based on the energy and frequency.
%   The function supports frequencies in both Hz and kHz.
%
% Example usage:
%   [E, cmax] = calcEnergyMaxCoeff(30, 5, 'Hz');   % for Hz
%   [E, cmax] = calcEnergyMaxCoeff(30, 5, 'kHz');  % for kHz
%

    % Convert frequency to Hz if it's in kHz
    if strcmp(unit, 'kHz')
        f = f * 1e3;  % Convert kHz to Hz
        tds = 1e-9;   % Time resolution in nanoseconds for kHz range
        timeRange = 50e-6;  % Time range in microseconds
    elseif strcmp(unit, 'Hz')
        tds = 1e-6;   % Time resolution in microseconds for Hz range
        timeRange = 50;   % Time range in seconds
    else
        error('Invalid unit. Use "Hz" or "kHz".');
    end
    
    % Set sampling frequency and time step
    Fs = 2 * f;             % Sampling frequency (Nyquist rate)
    ts = 1 / Fs;            % Sampling interval
    
    % Define time support based on the unit
    td = (-timeRange:tds:timeRange).';  % Time vector
    
    % Initialize the signal
    maxSig = zeros(length(td), 1);  % Preallocate signal
    
    % Coefficients for sinc function (ones for simplicity)
    an = ones(2 * N + 1, 1);
    
    % Build the signal as a sum of weighted sinc functions
    for n = -N:N
        maxSig = maxSig + an(n + N + 1) * sinc((td - n * ts) / ts);
    end
    
    % Calculate the energy of the signal
    E = sum(abs(maxSig).^2) * tds;
    
    % Calculate the maximum coefficient based on the energy and frequency
    cmax = sqrt(E * 2 * f); 
    %cmax calc is based on the paper:. A. Papoulis, “Limits on bandlimited signals,” Proceedings of the IEEE, vol. 55, no. 10, pp. 16771686, 1967.
    
end
