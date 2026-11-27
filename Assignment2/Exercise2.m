clc; clear; close all;
rng(22); % Reproducibility

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part 1: 2-PAM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ts = 1;
Eb = 1; % unitary bit energy (from Parseval's theorem)

SNR_dB = linspace(0, 12, 100);
Eb_N0  = 10.^(SNR_dB / 10);
Pb_analytical = 0.5 * erfc(sqrt(Eb_N0));

figure;
semilogy(SNR_dB, Pb_analytical, 'b-', 'LineWidth', 1.8);
xlabel('$\frac{E_b}{N_0}$ (dB)', 'Interpreter','latex')
ylabel('Bit Error Probability $P_b(e)$', 'Interpreter','latex')
title('\textbf{2-PAM - analytical BER}', 'Interpreter','latex')
grid on
set(gca, fontsize = 18)
grid on;


% Monte Carlo simulation
rng(22); % Reset seed for reproducibility

SNR_dB_mc     = 0:1:8;
target_errors = 100;
block_size    = 1e5;
BER_mc        = zeros(size(SNR_dB_mc));

for j = 1:length(SNR_dB_mc)

    % Convert SNR from dB to linear scale
    Eb_N0_lin = 10^(SNR_dB_mc(j) / 10);
    
    N0    = Eb / Eb_N0_lin;
    sigma = sqrt(N0 / 2); % noise standard deviation
    
    N_bit_errors    = 0;
    N_transmissions = 0;
    
    while N_bit_errors < target_errors

        % generate a vector of random bits {0, 1}
        x = randi([0, 1], block_size, 1);
        
        % apply binary labeling
        xb = (2 * x - 1) * sqrt(Eb);
        
        % add AWGN with variance N0/2
        noise = sigma * randn(block_size, 1);
        xb_r  = xb + noise;
        
        % minimum distance detection (sign decision)
        x_r = (xb_r >= 0);
        
        % count bit errors and total bits simulated
        N_bit_errors    = N_bit_errors + sum(x_r ~= x);
        N_transmissions = N_transmissions + block_size;
    end
    
    % estimate BER
    BER_mc(j) = N_bit_errors / N_transmissions;

    fprintf('Eb/N0 = %d dB, BER = %.6e, Errors = %d\n', ...
        SNR_dB_mc(j), BER_mc(j), N_bit_errors);
end

% Overlay Monte Carlo results onto the analytical plot
figure;
semilogy(SNR_dB, Pb_analytical, 'b-', 'LineWidth', 1.8);
hold on
semilogy(SNR_dB_mc, BER_mc, 'ro', 'MarkerSize', 8, 'LineWidth', 1.5);
xlabel('$\frac{E_b}{N_0}$ (dB)', 'Interpreter', 'latex');
ylabel('Bit Error Probability', 'Interpreter', 'latex');
title('\textbf{2-PAM -- Analytical vs. Monte Carlo BER}', 'Interpreter', 'latex');
legend('Analytical $P_b$', 'Monte Carlo BER', 'Interpreter', 'latex', 'Location', 'southwest');
ylim([1e-6, 1]); grid on;
set(gca, 'FontSize', 18);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part 2: DSSS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = 12; % spreading factor
Eb = N; % energy per bit
Ts = 1;

% generate random binary spreading sequence
rng(22); % reproducibility
c1 = randi([0 1], 1, N); % random binary spreading sequence

% convert binary sequence into bipolar sequence {-1,+1}
c1_bipolar = 2*c1 - 1;

disp('Binary spreading sequence:');
disp(c1);

disp('Bipolar spreading sequence:');
disp(c1_bipolar);


% SNR from 0 to 8 dB
SNR_dB_DSSS = 0:1:8;
Eb_N0_DSSS = 10.^(SNR_dB_DSSS/10); % linear value

% theoretical BER for 2-PAM
Pb_analytical_DSSS = 0.5 * erfc(sqrt(Eb_N0_DSSS));

figure;
semilogy(SNR_dB_DSSS, Pb_analytical_DSSS, ...
    'b-', 'LineWidth', 1.8);

xlabel('$E_b/N_0$ (dB)', 'Interpreter','latex');
ylabel('Bit Error Probability $P_b(e)$', 'Interpreter','latex');
title('\textbf{DSSS - Analytical BER}', 'Interpreter','latex');
grid on;
set(gca, 'FontSize', 18);


% Monte Carlo DSSS Simulation
rng(22); % reproducibility

target_errors = 100; % at least 100 errors per SNR
block_size = 1e5; % number of information bits per block

BER_DSSS = zeros(size(SNR_dB_DSSS));

for j = 1:length(SNR_dB_DSSS)

    % convert SNR from dB to linear scale
    Eb_N0_lin = 10^(SNR_dB_DSSS(j)/10);

    % calculate noise variance
    N0 = Eb / Eb_N0_lin;
    sigma = sqrt(N0/2);

    N_bit_errors = 0;
    N_transmissions = 0;

    while N_bit_errors < target_errors

        % generate random information bits
        x = randi([0 1], block_size, 1);

        % binary to bipolar mapping
        xb = 2*x - 1;

        % spread each information bit using the spreading sequence
        spread_signal = xb * c1_bipolar;

        % add independent Gaussian noise to each chip
        noise = sigma * randn(block_size, N);

        received_signal = spread_signal + noise;

        % despreading: correlate received chips with spreading sequence
        despread_signal = sum(received_signal .* c1_bipolar, 2); % project r onto c1

        % decision based on the sign of the correlation
        x_detected = (despread_signal >= 0);

        % count bit errors
        N_bit_errors = N_bit_errors + sum(x_detected ~= x);

        N_transmissions = N_transmissions + block_size;

    end

    % estimate BER
    BER_DSSS(j) = N_bit_errors / N_transmissions;

    fprintf('Eb/N0 = %d dB, BER = %.6e, Errors = %d\n', ...
        SNR_dB_DSSS(j), BER_DSSS(j), N_bit_errors);

end

% plot analytical vs simulated BER
figure;
semilogy(SNR_dB_DSSS, Pb_analytical_DSSS, 'b-', 'LineWidth', 1.8);
hold on;
semilogy(SNR_dB_DSSS, BER_DSSS, 'ro', 'MarkerSize', 8, 'LineWidth', 1.5);

xlabel('$E_b/N_0$ (dB)', 'Interpreter','latex');
ylabel('Bit Error Probability $P_b(e)$', 'Interpreter','latex');
title('\textbf{DSSS - Analytical vs Monte Carlo BER}', 'Interpreter','latex');
legend('Analytical BER', 'Monte Carlo BER', 'Location','southwest');
grid on;
ylim([1e-6 1]);
set(gca, 'FontSize', 18);
hold off