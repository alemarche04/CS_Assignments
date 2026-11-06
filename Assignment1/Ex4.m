clc
clear
close all

% =========================================================================
%% 4a: PERFECT BINARY (7,4) HAMMING CODE
% =========================================================================
fprintf('——— 4a: Perfect Hamming Code (7,4) ———\n');

% code Parameters
n1 = 7; k1 = 4; t1 = 1;

% generator matrix
I1 = eye(k1);
P1 = [1 1 0;
      1 0 1;
      0 1 1;
      1 1 1];
G1 = [I1 P1];

% generate all possible message vectors
V1 = dec2bin(0:2^k1-1) - '0';

% generate the codebook
C1 = mod(V1 * G1, 2);

% verify perfect code property (compare decision regions vs Hamming spheres)
Y1      = dec2bin(0:2^n1-1) - '0'; % all possible transmitted vectors
D_c1    = cell(size(C1, 1), 1);
S_c1    = cell(size(C1, 1), 1);

for j = 1:size(Y1, 1)
    y = Y1(j, :);
    d_H = sum(y ~= C1, 2);
    
    % minimum distance decision region
    d_min = min(d_H);
    for idx = find(d_H == d_min)'
        D_c1{idx} = [D_c1{idx}; y];
    end
    
    % Hamming sphere of radius t = 1
    for idx = find(d_H <= t1)'
        S_c1{idx} = [S_c1{idx}; y];
    end
end

% compare decision regions and Hamming spheres
perfect1 = true;
for i = 1:size(C1, 1)
    perfect1 = perfect1 && isequal(D_c1{i}, S_c1{i});
end

if perfect1
    fprintf('VERIFICATION: The (7,4) Hamming code is PERFECT.\n\n');
else
    fprintf('VERIFICATION: The (7,4) code is NOT perfect.\n\n');
end


% analytical curve of the Error Codeword Probability
% focus on small p (logarithmic scale)
p_small = logspace(-5, -0.5, 100);
P_c1_small = zeros(size(p_small));
for i = 0:t1
    P_c1_small = P_c1_small + nchoosek(n1, i) * p_small.^i .* (1-p_small).^(n1-i);
end
P_w1_small = 1 - P_c1_small;

figure;
loglog(p_small, P_w1_small, 'b-', 'LineWidth', 2);
grid on; axis tight;
xlabel('BSC crossover probability $p$', 'Interpreter', 'latex');
ylabel('$P_w(E)$', 'Interpreter', 'latex');
set(gca, 'FontSize', 18);


% Monte Carlo Simulation (small p interval to avoid excessive execution time)
p_vec = logspace(-2, -0.5, 12); % 12 points in [0.01, 0.1]
target_errors = 100;
alpha = 0.05; % parameter for confidence intervals

P_W_MC1     = zeros(1, length(p_vec));
CI_low1     = zeros(1, length(p_vec));
CI_high1    = zeros(1, length(p_vec));

rng(22); % reproducibility
for i = 1:length(p_vec)
    pc = p_vec(i);
    N_WE = 0; N_trans = 0;
    while N_WE < target_errors
        v = randi([0 1], 1, k1);
        c = mod(v * G1, 2);
        e = rand(1, n1) < pc;
        y = xor(c, e);
        
        d_H = sum(C1 ~= y, 2);
        [~, decoded_idx] = min(d_H);
        v_hat = V1(decoded_idx, :);
        
        N_WE = N_WE + any(v_hat ~= v);
        N_trans = N_trans + 1;
    end

    P_W_MC1(i)  = N_WE / N_trans;
    r = N_WE; f = N_trans - N_WE;
    CI_low1(i)  = betainv(alpha/2, r, f+1);
    CI_high1(i) = betainv(1-alpha/2, r+1, f);
end

% analytical curve for comparison
P_c1_th = zeros(size(p_vec));
for i = 0:t1
    P_c1_th = P_c1_th + nchoosek(n1, i) * p_vec.^i .* (1-p_vec).^(n1-i);
end
P_w1_th = 1 - P_c1_th; % exact for perfect code

% plot perfect code (7,4)
figure;
loglog(p_small, P_w1_small, 'b-', 'LineWidth', 2); hold on;
errorbar(p_vec, P_W_MC1, P_W_MC1 - CI_low1, CI_high1 - P_W_MC1, '-ro', ...
    'LineWidth', 1.5, 'MarkerSize', 8, 'MarkerFaceColor', 'r');
grid on; axis tight;
xlabel('BSC crossover probability $p$', 'Interpreter', 'latex');
ylabel('$P_w(E)$', 'Interpreter', 'latex');
legend('Analytical $P_w(E)$', 'Monte Carlo (100 word errors)', 'Interpreter', 'latex', 'Location', 'northwest');
set(gca, 'FontSize', 18);


% =========================================================================
%% 4c: NON-PERFECT BINARY (8,4) CODE
% =========================================================================
fprintf('——— 4c: Non-Perfect Binary Code (8,4) ———\n');

% code Parameters
n2 = 8; k2 = 4; t2 = 1;

% generator matrix
I2 = eye(k2);
P2 = [1 1 0 1;
      1 0 1 1;
      0 1 1 1;
      1 1 1 0];
G2 = [I2 P2];

% generate all possible information vectors
V2 = dec2bin(0:2^k2-1) - '0';

% generate the codebook
C2 = mod(V2 * G2, 2);

% verify non-perfect code property (compare decision regions vs Hamming spheres)
Y2      = dec2bin(0:2^n2-1) - '0';
D_c2    = cell(size(C2, 1), 1);
S_c2    = cell(size(C2, 1), 1);

for j = 1:size(Y2, 1)
    y = Y2(j, :); % all possible transmitted vectors
    d_H = sum(y ~= C2, 2);
    
    d_min = min(d_H);
    for idx = find(d_H == d_min)'
        D_c2{idx} = [D_c2{idx}; y];
    end
    
    for idx = find(d_H <= t2)'
        S_c2{idx} = [S_c2{idx}; y];
    end
end

perfect2 = true;
for i = 1:size(C2, 1)
    perfect2 = perfect2 && isequal(D_c2{i}, S_c2{i});
end

if perfect2
    fprintf('VERIFICATION: The (8,4) code is PERFECT.\n\n');
else
    fprintf('VERIFICATION: The (8,4) code is NOT perfect.\n\n');
end


% analytical curve of the Error Codeword Probability
% upper bound formula: 1 - sum_{i=0}^t nchoosek(n,i)*p^i*(1-p)^(n-i))
% construct Decision Region D(0) for exact performance calculation
E2      = dec2bin(0:2^n2-1, n2) - '0';
w_e2    = sum(E2, 2);
in_D02  = false(2^n2, 1);
for m = 1:2^n2
    d = sum(C2 ~= E2(m,:), 2);
    in_D02(m) = (w_e2(m) == min(d));
end
D02_idx = find(in_D02);

P_c2_small = zeros(size(p_small));
for i = 1:numel(D02_idx)
    m = D02_idx(i);
    P_c2_small = P_c2_small + p_small.^w_e2(m) .* (1-p_small).^(n2-w_e2(m));
end
P_w2_small = 1 - P_c2_small;

P_c2_bound_small = zeros(size(p_small));
for i = 0:t2
    P_c2_bound_small = P_c2_bound_small + nchoosek(n2, i) * p_small.^i .* (1-p_small).^(n2-i);
end
P_w2_bound_small = 1 - P_c2_bound_small; % Upper Bound formula for non-perfect code

figure;
loglog(p_small, P_w2_small, 'b-', 'LineWidth', 2);
grid on; axis tight;
xlabel('BSC crossover probability $p$', 'Interpreter', 'latex');
ylabel('$P_w(E)$', 'Interpreter', 'latex');
set(gca, 'FontSize', 18);


% Monte Carlo Simulation for C(8,4)
P_W_MC2     = zeros(1, length(p_vec));
CI_low2     = zeros(1, length(p_vec));
CI_high2    = zeros(1, length(p_vec));

rng(22); % reproducibility
for i = 1:length(p_vec)
    pc = p_vec(i);
    N_WE = 0; N_trans = 0;
    while N_WE < target_errors
        v = randi([0 1], 1, k2);
        c = mod(v * G2, 2);
        e = rand(1, n2) < pc;
        y = xor(c, e);
        
        d_H = sum(C2 ~= y, 2);
        idx_tx = v * (2.^(k2-1:-1:0))' + 1;
        
        % word error only when another codeword is strictly closer
        N_WE = N_WE + (d_H(idx_tx) > min(d_H));
        N_trans = N_trans + 1;
    end
    P_W_MC2(i)  = N_WE / N_trans;
    r = N_WE; f = N_trans - N_WE;
    CI_low2(i)  = betainv(alpha/2,   r,   f+1);
    CI_high2(i) = betainv(1-alpha/2, r+1, f);
end

% analytical curve for comparison
P_c2_bound = zeros(size(p_vec));
for i = 0:t2
    P_c2_bound = P_c2_bound + nchoosek(n2, i) * p_vec.^i .* (1-p_vec).^(n2-i);
end
P_w2_bound = 1 - P_c2_bound; % Upper Bound formula for non-perfect code

P_c2_exact = zeros(size(p_vec));
for i = 1:numel(D02_idx)
    m = D02_idx(i);
    P_c2_exact = P_c2_exact + p_vec.^w_e2(m) .* (1-p_vec).^(n2-w_e2(m));
end
P_w2_exact = 1 - P_c2_exact;

% plot non-perfect code (8,4)
figure;
loglog(p_small, P_w2_bound_small, 'k--', 'LineWidth', 1.8); hold on;
loglog(p_small, P_w2_small, 'b-', 'LineWidth', 2);
errorbar(p_vec, P_W_MC2, P_W_MC2 - CI_low2, CI_high2 - P_W_MC2, '-ro', ...
    'LineWidth', 1.5, 'MarkerSize', 8, 'MarkerFaceColor', 'r');
grid on; axis tight;
xlabel('BSC crossover probability $p$', 'Interpreter', 'latex');
ylabel('$P_w(E)$', 'Interpreter', 'latex');
legend('Sphere Upper Bound', 'Exact Analytical $P_w(E)$', ...
       'Monte Carlo (100 word errors)', 'Interpreter', 'latex', 'Location', 'northwest');
set(gca, 'FontSize', 18);


% =========================================================================
%% PERFORMANCE COMPARISON
% =========================================================================
figure;
loglog(p_small, P_w1_small, 'b-', 'LineWidth', 2); hold on;
loglog(p_small, P_w2_small, 'r-', 'LineWidth', 2);
loglog(p_small, P_w2_bound_small, 'k--', 'LineWidth', 1.5);
grid on; axis tight;
xlabel('BSC crossover probability $p$', 'Interpreter', 'latex');
ylabel('$P_w(E)$', 'Interpreter', 'latex');
legend('$\mathcal{C}(7,4)$ Analytical', ...
       '$\mathcal{C}(8,4)$ Exact Analytical', ...
       '$\mathcal{C}(8,4)$ Sphere Upper Bound', ...
       'Interpreter', 'latex', 'Location', 'southeast');
set(gca, 'FontSize', 18);