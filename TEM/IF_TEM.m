
% encodes the signal x using Integrate and Fire Time Encoding 
%
% Inputs:
%   - x: Signal (vector) to be encoded
%   - dt: Signal resolution (scalar)
%   - bias:  Bias (scalar)
%   - th: Threshold (scalar)
%   - k: Integrator accelerator

% Outputs:
%   - tr: Trigger time indices (vector)
%   - intg_o: Integrator output at each time step 

% Author: Aseel Omar
function [tr,intg_o] = IF_TEM(x,dt,bias,th,k)


    % Preprocessing
    len = length(x); j = 1;

    integOut = 0;

     % Define the function to compute integrator output at each step
    integOutCalc = @(integOut,i) integOut + dt*(bias+x(i))/k;

     % Main loop to process signal
    for i=1:len
        
      % Compute integrator output  
      integOut = integOutCalc(integOut,i);
      intg_o(i) =integOut;
      % Threshold firing condition
      if integOut >= th
        tr(j) = i; % Record trigger time index

        % Reset integrator and increment trigger counter
        j = j + 1;
        integOut = integOut - th;

      end
    end
       % Handle case where no triggers were recorded
    if isempty(tr)
        warning('No triggers were recorded; signal did not reach the threshold.');
    end 
 %   figure;plot(intg_y);
end