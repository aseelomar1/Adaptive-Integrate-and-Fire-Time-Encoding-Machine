%% Parameter Initialization
f = 10;  % Hz, frequency of the input signal
tds = 1e-6;  % Time resolution (s)
TP = 0.7;  % Total time length (s)
N = 2;  % Number of sinc components
th_c = 0.0188;  % IF-TEM threshold (adjustable)
kappa = 0.24;  % IF-TEM integrator acceleration
beta = 0.1;  % MAP parameter for AIF-TEM
bmin = 0.1;  % Minimum bias for AIF-TEM
alpha1 = 0.98;  % Alpha1 for MAP block
alpha2 = 0.3;  % Alpha2 for MAP block
w = 5;  % Window size for bias adaptation
biasbits =4; %bias allocated for the bias vetor for AIF-TEM
verbose = true;  % Enable verbose/debug mode

%% add the path of used functions
addpath('../HelperFunc');
addpath('../TEM');
%% Generate the Input Signal
%t = 0:td:TP;  % Time vector
[x, t, cn] = genSinc(f, tds, TP, N);
c = max(abs(x));  % Maximum amplitude

% Energy and Maximum Amplitude Calculation
[E_max, cmax] = calcEnergyMaxCoeff(f, N,'Hz');
% cmax calcuation is based on the paper: . A. Papoulis, “Limits on bandlimited signals,” Proceedings of the IEEE, vol. 55, no. 10, pp. 16771686, 1967.
%% IF-TEM Sampling
bias_if = cmax + beta;  % IF-TEM bias
r_c = ((kappa * th_c * f * 2) / (bias_if - cmax));

if verbose && (r_c > 1)
    disp('IF-TEM recovery condition does not hold');
end

[trIndices_if, y_out] = IF_TEM(x, tds, bias_if, th_c, kappa);
tn_if = t(trIndices_if)';
OS_c =  numel(trIndices_if)/((max(t)-min(t))/(1/(2*f)));% calcuate Oversampling of IF-TEM

%% IF-TEM Signal Recovery
xr_if = recover_TEM(t', tn_if, 2 * pi * f, bias_if, kappa, th_c);
if size(x) ~= size(xr_if)
    xr_if = xr_if';
end

cut = 0.5 * 10^5;
mse_if = immse(x(cut:end-cut), xr_if(cut:end-cut));
mse_db_if = 10 * log10(mse_if);

if verbose
    disp(['The sampling MSE in dB for IF-TEM: ', num2str(mse_db_if), ' with oversampling: ', num2str(OS_c)]);
end

%% AIF-TEM Sampling
bias_aif = bias_if;  % AIF-TEM bias initialization
L= 2^biasbits;
Delta_b = (bias_aif-bmin)/(L-1);  % Step size for bias adaptation
r_a = ((kappa * th_c * f * 2) / (beta));

if verbose && (r_a > 1)
    disp('AIF-TEM recovery condition does not hold');
end

[trIndices_aif, bias_adaptive, y_out_aif, c_n_aif] = AIF_TEM(x, tds, bias_aif, th_c, kappa, ...
    alpha1, alpha2, beta, Delta_b, bmin, w, verbose);
tn_aif = t(trIndices_aif)';
bias_adaptive=bias_adaptive(2:end-1);
OS_a =  numel(trIndices_aif)/((max(t)-min(t))/(1/(2*f)));% calcuate Oversampling of AIF-TEM

%% AIF-TEM Signal Recovery
xr_aif = recover_TEM(t', tn_aif, 2 * pi * f, bias_adaptive, kappa, th_c);
if size(x) ~= size(xr_aif)
    xr_aif = xr_aif';
end

mse_aif = immse(x(cut:end-cut), xr_aif(cut:end-cut));
mse_db_aif = 10 * log10(mse_aif);

if verbose
    disp(['The sampling MSE in dB for AIF-TEM: ', num2str(mse_db_aif), ' with oversampling: ', num2str(OS_a)]);
end

%% Plotting - Improved with Subplots for Comparison
figure;

subplot(2, 1, 1);
plot(t, x, 'b', 'LineWidth', 2);
hold on;
stem(tn_if, x(trIndices_if), 'r', 'filled');
plot(t, xr_if, 'g', 'LineWidth', 1.5);
legend('Original Signal', 'Sampled Points (IF-TEM)', 'Recovered Signal (IF-TEM)');
xlabel('Time (s)');
ylabel('Amplitude');
title('IF-TEM Sampling and Recovery');
grid on;

subplot(2, 1, 2);
plot(t, x, 'b', 'LineWidth', 2);
hold on;
stem(tn_aif, x(trIndices_aif), 'r', 'filled');
plot(t, xr_aif, 'g', 'LineWidth', 1.5);
legend('Original Signal', 'Sampled Points (AIF-TEM)', 'Recovered Signal (AIF-TEM)');
xlabel('Time (s)');
ylabel('Amplitude');
title('AIF-TEM Sampling and Recovery');
grid on;
%%%%%%%%%%%%%%%%



