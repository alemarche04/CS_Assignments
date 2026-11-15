clc
clear
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 4a: Perfect binary Hamming code C(7,4)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('——— 4a: Perfect Hamming Code (7,4) ———\n');
% code Parameters
n1 = 7; k1 = 4;
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

% Compute minimum Hamming distance
dmin1 = n1;
for i = 1:size(C1,1)
    for j = i+1:size(C1,1)
        d = sum(C1(i,:) ~= C1(j,:));
        dmin1 = min(dmin1, d);
    end
end

% Error correction capability
t1 = floor((dmin1-1)/2);
fprintf('Minimum distance dmin = %d\n', dmin1);
fprintf('Error correction capability t = %d\n', t1);

sphere_size1 = 2^k1 * sum(arrayfun(@(i) nchoosek(n1,i), 0:t1));
perfect1 = (sphere_size1 == 2^n1);
if perfect1
    fprintf('VERIFICATION: The (7,4) code is PERFECT.\n\n');
else
    fprintf('VERIFICATION: The (7,4) code is NOT perfect.\n\n');
end

% % verify perfect code property (compare decision regions vs Hamming spheres)
% Y1      = dec2bin(0:2^n1-1) - '0'; % all possible transmitted vectors
% D_c1    = cell(size(C1, 1), 1);
% S_c1    = cell(size(C1, 1), 1);
% for j = 1:size(Y1, 1)
%     y = Y1(j, :);
%     d_H = sum(y ~= C1, 2);
% 
%     % minimum distance decision region
%     d_min = min(d_H);
%     for idx = find(d_H == d_min)'
%         D_c1{idx} = [D_c1{idx}; y];
%     end
% 
%     % Hamming sphere of radius t = 1
%     for idx = find(d_H <= t1)'
%         S_c1{idx} = [S_c1{idx}; y];
%     end
% end
% 
% % compare decision regions and Hamming spheres
% perfect1 = true;
% for i = 1:size(C1, 1)
%     perfect1 = perfect1 && isequal(D_c1{i}, S_c1{i});
% end
% if perfect1
%     fprintf('VERIFICATION: The (7,4) Hamming code is PERFECT.\n\n');
% else
%     fprintf('VERIFICATION: The (7,4) code is NOT perfect.\n\n');
% end

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



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 4c: Non-perfect binary linear code C(8,4)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('——— 4c: Non-Perfect Binary Code (8,4) ———\n');

% code Parameters
n2 = 8; k2 = 4;

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

% Compute minimum Hamming distance
dmin2 = n2;
for i = 1:size(C2,1)
    for j = i+1:size(C2,1)
        d = sum(C2(i,:) ~= C2(j,:));
        dmin2 = min(dmin2, d);
    end
end

% Error correction capability
t2 = floor((dmin2-1)/2);
fprintf('Minimum distance dmin = %d\n', dmin2);
fprintf('Error correction capability t = %d\n', t2);

sphere_size2 = 2^k2 * sum(arrayfun(@(i) nchoosek(n2,i), 0:t2));
perfect2 = (sphere_size2 == 2^n2);
if perfect2
    fprintf('VERIFICATION: The (8,4) code is PERFECT.\n\n');
else
    fprintf('VERIFICATION: The (8,4) code is NOT perfect.\n\n');
end

% Analytical error probability: build disjoint translation-invariant D(0)
Y2 = dec2bin(0:2^n2-1) - '0';
w_Y2 = sum(Y2, 2);
[~, order_idx] = sort(w_Y2, 'ascend'); % prioritize smaller weights

assigned = false(2^n2, 1);
in_D02   = false(2^n2, 1);
pow2_n   = (2.^(n2-1:-1:0))';

for idx_s = 1:size(Y2, 1)
    j = order_idx(idx_s);
    if assigned(j)
        continue;
    end
    y = Y2(j, :);
    d_H = sum(C2 ~= y, 2);
    % Zero codeword is the first row of C2
    if d_H(1) == min(d_H)
        in_D02(j) = true;
        % Mark all shifted versions (c XOR y) as assigned to their respective D(c)
        for c_idx = 1:size(C2, 1)
            y_shifted_id = xor(C2(c_idx, :), y) * pow2_n + 1;
            assigned(y_shifted_id) = true;
        end
    end
end

% Count error patterns outside D(0) by Hamming weight
N2 = zeros(1, n2+1);
for j = 1:size(Y2, 1)
    if ~in_D02(j)
        w = w_Y2(j);
        N2(w+1) = N2(w+1) + 1;
    end
end
disp('Error patterns by weight:');
disp(N2);

% Monte Carlo Simulation
P_W_MC2 = zeros(1,length(p_vec));
CI_low2 = zeros(1,length(p_vec));
CI_high2 = zeros(1,length(p_vec));
rng(22);
for i = 1:length(p_vec)
    pc = p_vec(i);
    N_WE = 0;
    N_trans = 0;
    while N_WE < target_errors
        % Random information vector
        idx_tx = randi(2^k2);
        v = V2(idx_tx,:);
        % Encoding
        c = mod(v * G2, 2);
        % BSC channel
        e = rand(1,n2) < pc;
        y = xor(c,e);
        % Minimum distance decoding (select a single codeword, breaking ties blindly)
        d_H = sum(C2 ~= y,2);
        d_min = min(d_H);
        nearest = find(d_H == d_min);
        decoded_idx = nearest(randi(numel(nearest)));
        if decoded_idx ~= idx_tx
            N_WE = N_WE + 1;
        end
        N_trans = N_trans + 1;
    end
    P_W_MC2(i) = N_WE/N_trans;
    % Confidence intervals
    r = N_WE;
    f = N_trans-N_WE;
    CI_low2(i) = betainv(alpha/2,r,f+1);
    CI_high2(i) = betainv(1-alpha/2,r+1,f);
end

% Analytical curve
P_W2_th = zeros(size(p_small));
for i = 0:n2
    P_W2_th = P_W2_th + N2(i+1) * p_small.^i .* (1-p_small).^(n2-i);
end

% Upper bound on the correct decoding probability
P_c2_bound = zeros(size(p_small));
for i = 0:t2
    P_c2_bound = P_c2_bound + ...
        nchoosek(n2, i) * p_small.^i .* (1 - p_small).^(n2 - i);
end
% Upper bound on the word error probability
P_w2_bound = 1 - P_c2_bound;

figure;
loglog(p_small, P_W2_th, 'b-', 'LineWidth', 2);
hold on;
loglog(p_small, P_w2_bound, 'k--', 'LineWidth', 2);
errorbar(p_vec, P_W_MC2, ...
    P_W_MC2-CI_low2, CI_high2-P_W_MC2, ...
    '-ro', 'LineWidth', 1.5, ...
    'MarkerSize', 8, 'MarkerFaceColor', 'r');
grid on;
axis tight;
xlabel('BSC crossover probability $p$', 'Interpreter', 'latex');
ylabel('$P_W(E)$', 'Interpreter', 'latex');
legend('Exact analytical $P_W(E)$', ...
       'Upper bound on $P_W(E)$', ...
       'Monte Carlo (100 word errors)', ...
       'Interpreter', 'latex', ...
       'Location', 'northwest');
set(gca, 'FontSize', 18);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Comparison between (7,4) and (8,4) codes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
loglog(p_small, P_w1_small, 'b-', 'LineWidth', 2);
hold on;
loglog(p_small, P_W2_th, 'r--', 'LineWidth', 2);
grid on;
axis tight;
xlabel('BSC crossover probability $p$', 'Interpreter', 'latex');
ylabel('$P_W(E)$', 'Interpreter', 'latex');
legend('Hamming (7,4)', ...
       'Extended Hamming (8,4)', ...
       'Interpreter', 'latex', ...
       'Location', 'northwest');
set(gca, 'FontSize', 18);