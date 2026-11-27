clear; close all; clc;

% Physical Constants
c   = 299792458;          % Speed of light [m/s]
k_B = 1.380649e-23;       % Boltzmann constant [J/K]
R_E = 6371e3;             % Earth radius [m]
h   = 35786e3;            % GEO satellite altitude [m]

% Selected System Parameters
f0      = 12e9;           % Carrier frequency [Hz] (12 GHz)
B       = 36e6;           % Channel bandwidth [Hz] (36 MHz)
P_TX    = 100;            % Satellite transmitted power [W] (20 dBW)
D_TX    = 1.0;            % Satellite TX antenna diameter [m]
eta_TX  = 0.65;           % Satellite TX antenna efficiency [-]
D_RX    = 0.6;            % User RX antenna diameter [m]
eta_RX  = 0.65;           % User RX antenna efficiency [-]
theta_e = 30;             % Elevation angle [deg]
T_sys   = 150;            % System noise temperature [K]
L_add   = 3.0;            % Additional losses [dB]
R_HD    = 5e6;            % Average bit rate per HD TV channel [bit/s]

lambda  = c / f0;         % Carrier wavelength [m]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part 1: Satellite-to-User distance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Question 2: Slant range at selected elevation angle (30 deg)
theta_sel = theta_e;      % [deg]
d_sel = sqrt((R_E + h)^2 - R_E^2 * cosd(theta_sel)^2) - R_E * sind(theta_sel);
fprintf('Part 1 - Slant range at theta_e = %.1f deg: %.2f km\n', theta_sel, d_sel / 1e3);

% Question 3: Slant range and FSPL across 10 deg <= theta_e <= 90 deg
theta_range = 10:1:90;
theta_rad   = deg2rad(theta_range);     % Elevation angle vector [rad]
d    = sqrt((R_E + h)^2 - R_E^2 * cos(theta_rad).^2) - R_E * sin(theta_rad); % [m]
L_FS = 20 * log10(4 * pi * d / lambda);                                      % [dB]

figure('Name', 'Part 1: Slant Range and FSPL', 'Color', 'w');
yyaxis left;
plot(theta_range, d / 1e3, 'LineWidth', 1.8);
ylabel('Slant Range d [km]');
grid on;

yyaxis right;
plot(theta_range, L_FS, '--', 'LineWidth', 1.8);
ylabel('Free-Space Path Loss L_{FS} [dB]');

xlabel('Elevation Angle \theta_e [deg]');
title('Slant Range and Free-Space Path Loss vs. Elevation Angle');
legend('Distance d [km]', 'L_{FS} [dB]', 'Location', 'northeast');
set(gca, 'FontSize', 16);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part 2: Link budget
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Question 5: TX and RX antenna gains [dBi]
G_TX_lin = eta_TX * (pi * D_TX / lambda)^2;
G_RX_lin = eta_RX * (pi * D_RX / lambda)^2;
G_TX_dBi = 10 * log10(G_TX_lin);
G_RX_dBi = 10 * log10(G_RX_lin);

% Question 6: EIRP [dBW]
P_TX_dBW = 10 * log10(P_TX);
EIRP_dBW = P_TX_dBW + G_TX_dBi;

% Questions 7 & 8: Free-space path loss and total losses [dB]
L_FS_dB  = 20 * log10(4 * pi * d_sel / lambda);
L_tot_dB = L_FS_dB + L_add;

% Question 9: Received power [dBW]
P_RX_dBW = EIRP_dBW - L_tot_dB + G_RX_dBi;

% Question 10: Thermal noise power [dBW]
P_N_lin  = k_B * T_sys * B;
P_N_dBW  = 10 * log10(P_N_lin);

% Question 11: Signal-to-Noise Ratio [dB]
SNR_dB   = P_RX_dBW - P_N_dBW;
SNR_lin  = 10^(SNR_dB / 10);

% Question 12: Print Link Budget Summary Table
fprintf('\n====================================================================\n');
fprintf('                  PART 2: LINK BUDGET SUMMARY TABLE                 \n');
fprintf('====================================================================\n');
fprintf('Carrier Frequency (f0)           : %10.2f GHz\n', f0 / 1e9);
fprintf('Bandwidth (B)                    : %10.2f MHz\n', B / 1e6);
fprintf('Elevation Angle (theta_e)        : %10.2f deg\n', theta_e);
fprintf('Slant Range (d)                  : %10.2f km\n', d_sel / 1e3);
fprintf('--------------------------------------------------------------------\n');
fprintf('Satellite TX Power (P_TX)        : %+10.2f dBW  (%.1f W)\n', P_TX_dBW, P_TX);
fprintf('Satellite TX Antenna Gain (G_TX) : %+10.2f dBi\n', G_TX_dBi);
fprintf('EIRP                             : %+10.2f dBW\n', EIRP_dBW);
fprintf('Free-Space Path Loss (L_FS)      : %10.2f dB\n', L_FS_dB);
fprintf('Additional Losses (L_add)        : %10.2f dB\n', L_add);
fprintf('Total Losses (L_tot)             : %10.2f dB\n', L_tot_dB);
fprintf('User RX Antenna Gain (G_RX)      : %+10.2f dBi\n', G_RX_dBi);
fprintf('Received Power (P_RX)            : %+10.2f dBW\n', P_RX_dBW);
fprintf('Noise Power (P_N)                : %+10.2f dBW\n', P_N_dBW);
fprintf('Signal-to-Noise Ratio (SNR)      : %+10.2f dB   (%.2f linear)\n', SNR_dB, SNR_lin);
fprintf('====================================================================\n\n');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part 3: Bit Rate & Shannon Capacity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Questions 13 & 14: Select highest spectral efficiency MODCOD (alpha = 0)
[eta_sel, modcod_sel, req_EsN0_sel] = select_best_modcod(SNR_dB);
Rs = B;                     % Symbol rate [symbol/s] (alpha = 0 => Rs = B)
Rb = Rs * eta_sel;          % Achievable bit rate [bit/s]

% Question 15: Simultaneous 5 Mbit/s HD TV channels
N_HD = floor(Rb / R_HD);

% Questions 16 & 17: Shannon Capacity and comparison
C_shannon = B * log2(1 + SNR_lin);

fprintf('====================================================================\n');
fprintf('           PART 3: BIT RATE & SHANNON CAPACITY COMPARISON           \n');
fprintf('====================================================================\n');
fprintf('Available Es/N0 (= SNR)          : %10.2f dB\n', SNR_dB);
fprintf('Selected DVB-S2X MODCOD          : %s (Req. Es/N0 = %.2f dB)\n', modcod_sel, req_EsN0_sel);
fprintf('Clear-Sky Link Margin            : %10.2f dB\n', SNR_dB - req_EsN0_sel);
fprintf('Spectral Efficiency (eta)        : %10.6f bit/symbol\n', eta_sel);
fprintf('Achievable Bit Rate (Rb)         : %10.2f Mbit/s\n', Rb / 1e6);
fprintf('Simultaneous HD Channels (5 Mbps): %10d channels\n', N_HD);
fprintf('Shannon Capacity (C)             : %10.2f Mbit/s\n', C_shannon / 1e6);
fprintf('Efficiency vs. Shannon (Rb / C)  : %10.2f %%\n', 100 * Rb / C_shannon);
fprintf('====================================================================\n\n');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part 4: System Trade-offs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Question 19: Bit rate vs. Transmitted Power P_TX (B = 36 MHz fixed)
P_TX_sweep  = linspace(5, 250, 1500);   % 5 W to 250 W
Rb_vs_P     = zeros(size(P_TX_sweep));
C_vs_P      = zeros(size(P_TX_sweep));

for i = 1:length(P_TX_sweep)
    eirp_i     = 10*log10(P_TX_sweep(i)) + G_TX_dBi;
    prx_i      = eirp_i - L_tot_dB + G_RX_dBi;
    snr_db_i   = prx_i - P_N_dBW;
    eta_i      = select_best_modcod(snr_db_i);
    Rb_vs_P(i) = B * eta_i;
    C_vs_P(i)  = B * log2(1 + 10^(snr_db_i / 10));
end

figure('Name', 'Part 4a: Bit Rate vs TX Power', 'Color', 'w');
plot(P_TX_sweep, Rb_vs_P / 1e6, 'b-', 'LineWidth', 2); hold on;
plot(P_TX_sweep, C_vs_P / 1e6, 'k--', 'LineWidth', 1.5);
plot(P_TX, Rb / 1e6, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
grid on;
xlabel('Satellite Transmitted Power P_{TX} [W]', 'FontWeight', 'bold');
ylabel('Bit Rate [Mbit/s]', 'FontWeight', 'bold');
title('DVB-S2X Bit Rate & Shannon Capacity vs. TX Power (B = 36 MHz)');
legend('DVB-S2X Bit Rate R_b', 'Shannon Capacity C', 'Baseline Operating Point (100 W)', ...
    'Location', 'southeast');
set(gca, 'FontSize', 16);

% Question 20: Bit rate vs. Bandwidth B (P_TX = 100 W fixed)
B_sweep = linspace(5e6, 100e6, 1500);  % 5 MHz to 100 MHz
Rb_vs_B = zeros(size(B_sweep));
C_vs_B  = zeros(size(B_sweep));

for i = 1:length(B_sweep)
    pn_dbw_i   = 10 * log10(k_B * T_sys * B_sweep(i));
    snr_db_i   = P_RX_dBW - pn_dbw_i;
    eta_i      = select_best_modcod(snr_db_i);
    Rb_vs_B(i) = B_sweep(i) * eta_i;
    C_vs_B(i)  = B_sweep(i) * log2(1 + 10^(snr_db_i / 10));
end

figure('Name', 'Part 4b: Bit Rate vs Bandwidth', 'Color', 'w');
plot(B_sweep / 1e6, Rb_vs_B / 1e6, 'b-', 'LineWidth', 2); hold on;
plot(B_sweep / 1e6, C_vs_B / 1e6, 'k--', 'LineWidth', 1.5);
plot(B / 1e6, Rb / 1e6, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
grid on;
xlabel('Channel Bandwidth B [MHz]', 'FontWeight', 'bold');
ylabel('Bit Rate [Mbit/s]', 'FontWeight', 'bold');
title('DVB-S2X Bit Rate & Shannon Capacity vs. Bandwidth (P_{TX} = 100 W)');
legend('DVB-S2X Bit Rate R_b', 'Shannon Capacity C', 'Baseline Operating Point (36 MHz)', ...
    'Location', 'northwest');
set(gca, 'FontSize', 16);


%% Optimal MODCOD Selection (ETSI EN 302 307-2)
function [best_eta, modcod, best_req_EsN0] = select_best_modcod(avail_EsN0)
    
    modcod_names = { ...
        'QPSK 2/9',        'QPSK 13/45',      'QPSK 9/20',       'QPSK 11/20', ...
        '8APSK 5/9-L',     '8APSK 26/45-L',   '8PSK 23/36',      '8PSK 25/36', ...
        '8PSK 13/18',      '16APSK 1/2-L',    '16APSK 8/15-L',   '16APSK 5/9-L', ...
        '16APSK 26/45',    '16APSK 3/5',      '16APSK 3/5-L',    '16APSK 28/45', ...
        '16APSK 23/36',    '16APSK 2/3-L',    '16APSK 25/36',    '16APSK 13/18', ...
        '16APSK 7/9',      '16APSK 77/90',    '32APSK 2/3-L',    '32APSK 32/45', ...
        '32APSK 11/15',    '32APSK 7/9',      '64APSK 32/45-L',  '64APSK 11/15', ...
        '64APSK 7/9',      '64APSK 4/5',      '64APSK 5/6',      '128APSK 3/4', ...
        '128APSK 7/9',     '256APSK 29/45-L', '256APSK 2/3-L',   '256APSK 31/45-L', ...
        '256APSK 32/45',   '256APSK 11/15-L', '256APSK 3/4'};
    
    spec_eff = [ ...
        0.434841, 0.567805, 0.889135, 1.088581, ...
        1.647211, 1.713601, 1.896173, 2.062148, ...
        2.145136, 1.972253, 2.104850, 2.193247, ...
        2.281645, 2.370043, 2.370043, 2.458441, ...
        2.524739, 2.635236, 2.745734, 2.856231, ...
        3.077225, 3.386618, 3.291954, 3.510192, ...
        3.620536, 3.841226, 4.206428, 4.338659, ...
        4.603122, 4.735354, 4.936639, 5.163248, ...
        5.355556, 5.065690, 5.241514, 5.417338, ...
        5.593162, 5.768987, 5.900855];
    
    req_EsN0 = [ ...
        -2.85, -2.03,  0.22,  1.45, ...
         4.73,  5.13,  6.12,  7.02, ...
         7.49,  5.97,  6.55,  6.84, ...
         7.51,  7.80,  7.41,  8.10, ...
         8.38,  8.43,  9.27,  9.71, ...
        10.65, 11.99, 11.10, 11.75, ...
        12.17, 13.05, 13.98, 14.81, ...
        15.47, 15.87, 16.55, 17.73, ...
        18.53, 16.98, 17.24, 18.10, ...
        18.59, 18.84, 19.57];
    
    feasible_indices = find(req_EsN0 <= avail_EsN0);
    if isempty(feasible_indices)
        best_eta      = 0;
        modcod        = 'None (Outage)';
        best_req_EsN0 = NaN;
    else
        [best_eta, rel_idx] = max(spec_eff(feasible_indices));
        best_idx            = feasible_indices(rel_idx);
        modcod              = modcod_names{best_idx};
        best_req_EsN0       = req_EsN0(best_idx);
    end
end