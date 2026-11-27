clc; clear; close all;
rng(22); % reproducibility

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part 1: 2-PAM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ts = 1;
Eb = 1; % unitary bit energy

SNR_dB   = linspace(0, 12, 200);
Eb_N0    = 10.^(SNR_dB / 10);
Pb_analytical = 0.5 * erfc(sqrt(Eb_N0));

% plot analytical BER (0 to 12 dB)
figure;
semilogy(SNR_dB, Pb_analytical, 'b-', 'LineWidth', 1.8);
xlabel('$E_b/N_0$ (dB)', 'Interpreter', 'latex');
ylabel('Bit Error Probability $P_b(e)$', 'Interpreter', 'latex');
title('\textbf{2-PAM - Analytical BER}', 'Interpreter', 'latex');
grid on;
ylim([1e-8 1]); grid on; xlim([0 12]); 
set(gca, 'FontSize', 16);

% Monte Carlo simulation (0 to 8 dB, step 1 dB, >= 100 errors)
rng(22);
SNR_dB_MC     = 0:1:8;
target_errors = 100;
block_size    = 1e5;
BER_MC        = zeros(size(SNR_dB_MC));

fprintf('\n___ Part 1: 2-PAM Monte Carlo ___\n');
for j = 1:length(SNR_dB_MC)
    Eb_N0_lin = 10^(SNR_dB_MC(j) / 10); % from db to linear scale
    N0        = Eb / Eb_N0_lin; % 
    sigma     = sqrt(N0 / 2); % noise standard deviation

    N_bit_errors    = 0;
    N_transmissions = 0;

    while N_bit_errors < target_errors
        x     = randi([0, 1], block_size, 1); % random bits (transmission)
        xb    = (2 * x - 1) * sqrt(Eb); % binary labeling
        noise = sigma * randn(block_size, 1); % independent noise

        xb_r  = xb + noise; % received signal
        x_r   = (xb_r >= 0); % minimum distance detection (sign detection)

        N_bit_errors    = N_bit_errors + sum(x_r ~= x);
        N_transmissions = N_transmissions + block_size;
    end

    BER_MC(j) = N_bit_errors / N_transmissions;
    fprintf('Eb/N0 = %d dB, BER = %.6e\n', SNR_dB_MC(j), BER_MC(j));
end

% overlay Monte Carlo onto analytical curve
figure;
semilogy(SNR_dB, Pb_analytical, 'b-', 'LineWidth', 1.8);
hold on;
semilogy(SNR_dB_MC, BER_MC, 'ro', 'MarkerSize', 8, 'LineWidth', 1.5);
xlabel('$E_b/N_0$ (dB)', 'Interpreter', 'latex');
ylabel('Bit Error Probability $P_b(e)$', 'Interpreter', 'latex');
title('\textbf{2-PAM - Analytical vs. Monte Carlo BER}', 'Interpreter', 'latex');
legend('Analytical $P_b$', 'Monte Carlo BER', 'Interpreter', 'latex', 'Location', 'southwest');
ylim([1e-8 1]); grid on; xlim([0 12]); grid on;
set(gca, 'FontSize', 18);
hold off;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part 2: Single-User DSSS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N  = 12; % spreading factor
Eb = N;  % energy per bit (Tc = 1)
Tc = 1;

% generate random binary spreading sequence
rng(22);
c1         = randi([0 1], 1, N); % random binary sequence of length N
c1_bipolar = 2*c1 - 1; % bipolar binary sequence

SNR_dB              = linspace(0, 12, 200);
Eb_N0               = 10.^(SNR_dB / 10);
SNR_dB_DSSS         = 0:1:8;
Eb_N0_DSSS          = 10.^(SNR_dB_DSSS / 10); % linear SNR value
Pb_analytical_DSSS  = 0.5 * erfc(sqrt(Eb_N0)); % analytical error probability

% plot analytical DSSS curve
figure;
semilogy(SNR_dB, Pb_analytical_DSSS, 'b-', 'LineWidth', 1.8);
xlabel('$E_b/N_0$ (dB)', 'Interpreter', 'latex');
ylabel('Bit Error Probability $P_b(e)$', 'Interpreter', 'latex');
title('\textbf{DSSS - Analytical BER}', 'Interpreter', 'latex');
grid on;
ylim([1e-8 1]); grid on; xlim([0 12]); 
set(gca, 'FontSize', 18);

% Monte Carlo DSSS Simulation
rng(22);
BER_DSSS = zeros(size(SNR_dB_DSSS));

fprintf('\n\n___ Part 2: Single-User DSSS Monte Carlo ___\n');
for j = 1:length(SNR_dB_DSSS)
    Eb_N0_lin = 10^(SNR_dB_DSSS(j) / 10);
    N0        = Eb / Eb_N0_lin;
    sigma     = sqrt(N0 / 2);

    N_bit_errors    = 0;
    N_transmissions = 0;

    while N_bit_errors < target_errors
        x               = randi([0 1], block_size, 1); % vector of random bits
        xb              = 2*x - 1; % binary labeling
        spread_signal   = xb * c1_bipolar; 
        noise           = sigma * randn(block_size, N);
        received_signal = spread_signal + noise;

        despread_signal = sum(received_signal .* c1_bipolar, 2); % project r onto c1
        x_detected      = (despread_signal >= 0);

        N_bit_errors    = N_bit_errors + sum(x_detected ~= x);
        N_transmissions = N_transmissions + block_size;
    end

    BER_DSSS(j) = N_bit_errors / N_transmissions;
    fprintf('Eb/N0 = %d dB, BER = %.6e\n', SNR_dB_DSSS(j), BER_DSSS(j));
end

% overlay simulated curve onto analytical curve
figure;
semilogy(SNR_dB, Pb_analytical_DSSS, 'b-', 'LineWidth', 1.8);
hold on;
semilogy(SNR_dB_DSSS, BER_DSSS, 'ro', 'MarkerSize', 8, 'LineWidth', 1.5);
xlabel('$E_b/N_0$ (dB)', 'Interpreter', 'latex');
ylabel('Bit Error Probability $P_b(e)$', 'Interpreter', 'latex');
title('\textbf{DSSS - Analytical vs. Monte Carlo BER}', 'Interpreter', 'latex');
legend('Analytical BER', 'Monte Carlo BER', 'Interpreter', 'latex', 'Location', 'southwest');
ylim([1e-8 1]); grid on; xlim([0 12]); 
set(gca, 'FontSize', 18);
hold off;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part 3: CDMA Downlink
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N  = 12;
Eb = N;

p_values    = [0, 2, 4, 12];
SNR_dB_sim  = 0:1:8;
SNR_dB_plot = linspace(0, 20, 400);
Eb_N0_plot  = 10.^(SNR_dB_plot / 10);
Pb_2PAM_ref = 0.5 * erfc(sqrt(Eb_N0_plot));

% % plot 2-PAM curve
% figure;
% semilogy(SNR_dB_plot, Pb_2PAM_ref, 'r--', 'LineWidth', 1.5); 
% xlabel('$E_b/N_0$ [dB]', 'Interpreter', 'latex');
% ylabel('BER', 'Interpreter', 'latex');
% title(sprintf('\\textbf{CDMA downlink}'), 'Interpreter', 'latex');
% legend('2-PAM', 'Location', 'northeast');
% grid on;
% ylim([1e-8 1]); xlim([0, 12]);
% set(gca,'FontSize',18);

% fixed sequence pairs for Downlink: p = |<c1', c2'>|
c1_DL_bipolar = ones(1, N);
c2_DL_list = [
    1  1  1  1  1  1 -1 -1 -1 -1 -1 -1;  % p = 0
    1  1  1  1  1  1  1 -1 -1 -1 -1 -1;  % p = 2
    1  1  1  1  1  1  1  1 -1 -1 -1 -1;  % p = 4
    1  1  1  1  1  1  1  1  1  1  1  1   % p = 12
    ];

fprintf('\n\n___ Part 3: CDMA Downlink ___\n');
for idx = 1:length(p_values)
    p        = p_values(idx);
    c1_b     = c1_DL_bipolar;
    c2_b     = c2_DL_list(idx, :);
    rho_12   = sum(c1_b .* c2_b); % cross correlation btw c1_b and c2_b

    % exact analytical curve
    Pb_exact_DL = 0.25 * erfc((1 + p/N) * sqrt(Eb_N0_plot)) + ...
        0.25 * erfc((1 - p/N) * sqrt(Eb_N0_plot));

    % approximate curve
    if p == 0
        Pb_approx_DL = 0.5 * erfc(sqrt(Eb_N0_plot));
    else
        Pb_approx_DL = 0.25 * erfc((1 - p/N) * sqrt(Eb_N0_plot));
    end

    % Monte Carlo simulation (0 to 8 dB)
    rng(22);
    BER_sim_DL = zeros(size(SNR_dB_sim));

    for j = 1:length(SNR_dB_sim)
        Eb_N0_lin = 10^(SNR_dB_sim(j) / 10);
        N0        = Eb / Eb_N0_lin;
        sigma     = sqrt(N0 / 2); % noise standard deviation

        N_bit_errors    = 0;
        N_transmissions = 0;

        while N_bit_errors < target_errors
            x1  = randi([0 1], block_size, 1);
            x2  = randi([0 1], block_size, 1);
            xb1 = 2*x1 - 1;
            xb2 = 2*x2 - 1;

            received_signal = xb1 * c1_b + xb2 * c2_b + sigma * randn(block_size, N);
            despread_signal = sum(received_signal .* c1_b, 2);
            x_detected      = (despread_signal >= 0);

            % compare with x1 (User 1's bits)
            N_bit_errors    = N_bit_errors + sum(x_detected ~= x1);
            N_transmissions = N_transmissions + block_size;
        end

        BER_sim_DL(j) = N_bit_errors / N_transmissions;
        fprintf('Eb/N0 = %d dB, BER = %.6e\n', SNR_dB_sim(j), BER_sim_DL(j));
    end

    % plot 2-PAM, Exact Analytic, Simulation, and Approximated curves
    figure;
    semilogy(SNR_dB_plot, Pb_2PAM_ref, 'r--', 'LineWidth', 1.5); 
    hold on;
    semilogy(SNR_dB_plot, Pb_exact_DL, 'b-', 'LineWidth', 1.5);
    semilogy(SNR_dB_sim,  BER_sim_DL, 'ko-', 'MarkerSize', 7, 'LineWidth', 1.5);
    semilogy(SNR_dB_plot, Pb_approx_DL, 'm--', 'LineWidth', 1.5);
    xlabel('$E_b/N_0$ [dB]', 'Interpreter', 'latex');
    ylabel('BER', 'Interpreter', 'latex');
    title(sprintf('\\textbf{CDMA downlink (N = %d, p = %d)}', N, p), 'Interpreter', 'latex');
    legend('2-PAM', 'analytic', 'simulation', 'approximated', 'Location', 'northeast');
    
    grid on;
    ylim([1e-8 1]); xlim([0, 12]);
    set(gca,'FontSize',18);
    hold off;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part 4: CDMA Uplink
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% fixed sequence pairs for Uplink: p = max_{m in 0..N-1} |rho(m)|
c1_UL_bipolar = [ 1  1  1 -1  1 -1 -1 -1 -1  1 -1  1];
c2_UL_list = [
     1  1 -1  1 -1 -1  1  1 -1  1 -1 -1;   % p = 0,  M_p = 12 (rho(m)=0 for all m)
     1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1;   % p = 2,  M_p = 12
     1  1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1;   % p = 4,  M_p = 6
     1  1  1 -1  1 -1 -1 -1 -1  1 -1  1 ]; % p = 12, M_p = 2  (c2 = c1)


fprintf('\n\n___ Part 4: CDMA Uplink (Cyclic Shifts) ___\n');
for idx = 1:length(p_values)
    p    = p_values(idx); % peak cross-correlation
    c1_b = c1_UL_bipolar;
    c2_b = c2_UL_list(idx, :);

    % compute periodic cross-correlation for all shifts m = 0,...,N-1
    rho_shifts = zeros(1, N);
    for m = 0:N-1
        c2_shifted        = circshift(c2_b, [0, m]);
        rho_shifts(m + 1) = sum(c1_b .* c2_shifted);
    end

    M_p     = sum(abs(rho_shifts) == p); % # of shifts for which rho(m) = p

    fprintf('\n[Uplink p = %d, M_p = %d]', p, M_p);

    % exact analytical Uplink BER averaged over all N cyclic shifts
    Pb_exact_UL = zeros(size(Eb_N0_plot));
    for m = 1:N
        rho_m = rho_shifts(m);
        Pb_exact_UL = Pb_exact_UL + (1 / N) * ( ...
            0.25 * erfc((1 + rho_m/N) * sqrt(Eb_N0_plot)) + ...
            0.25 * erfc((1 - rho_m/N) * sqrt(Eb_N0_plot)) );
    end

    % approximate Uplink BER based on p and multiplicity M_p
    if p == 0
        Pb_approx_UL = 0.5 * erfc(sqrt(Eb_N0_plot));
    else
        Pb_approx_UL = (M_p / (4 * N)) * erfc((1 - p/N) * sqrt(Eb_N0_plot));
    end

    % Monte Carlo simulation with independent uniform cyclic shift per bit
    rng(22);
    BER_sim_UL = zeros(size(SNR_dB_sim));

    % precompute all N cyclic shifts of c2_b as rows of an (N x N) matrix
    C2_shifts = zeros(N, N);
    for m = 0:N-1
        C2_shifts(m + 1, :) = circshift(c2_b, [0, m]);
    end

    for j = 1:length(SNR_dB_sim)
        Eb_N0_lin = 10^(SNR_dB_sim(j) / 10);
        N0        = Eb / Eb_N0_lin;
        sigma     = sqrt(N0 / 2);

        N_bit_errors    = 0;
        N_transmissions = 0;

        while N_bit_errors < target_errors
            x1  = randi([0 1], block_size, 1);
            x2  = randi([0 1], block_size, 1);
            xb1 = 2*x1 - 1;
            xb2 = 2*x2 - 1;

            % independent random cyclic shift m in {0, ..., N-1} for each bit
            shift_indices   = randi([1, N], block_size, 1);
            spread_signal2  = xb2 .* C2_shifts(shift_indices, :);

            received_signal = xb1 * c1_b + spread_signal2 + sigma * randn(block_size, N);
            despread_signal = sum(received_signal .* c1_b, 2);
            x_detected      = (despread_signal >= 0);

            N_bit_errors    = N_bit_errors + sum(x_detected ~= x1);
            N_transmissions = N_transmissions + block_size;
        end

        BER_sim_UL(j) = N_bit_errors / N_transmissions;
        fprintf('Eb/N0 = %d dB, BER = %.6e\n', SNR_dB_sim(j), BER_sim_UL(j));
    end

    % plot 2-PAM, Exact Analytic, Simulation, and Approximated curves
    figure;
    semilogy(SNR_dB_plot, Pb_2PAM_ref, 'r--', 'LineWidth', 1.5); hold on;
    semilogy(SNR_dB_plot, Pb_exact_UL, 'b-',  'LineWidth', 1.5);
    semilogy(SNR_dB_sim, BER_sim_UL, 'ko-', 'MarkerSize', 7, 'LineWidth', 1.5);
    semilogy(SNR_dB_plot, Pb_approx_UL, 'm--', 'LineWidth', 1.5);
    xlabel('$E_b/N_0$ [dB]', 'Interpreter', 'latex');
    ylabel('BER', 'Interpreter', 'latex');
    title(sprintf('\\textbf{CDMA uplink (N = %d, p = %d)}', N, p), 'Interpreter', 'latex');
    legend('2-PAM', 'analytic', 'simulation', 'approximated', 'Location', 'northeast');
    ylim([1e-8 1]); xlim([0, 12]);
    grid on;
    set(gca, 'FontSize', 18);
    hold off;
end