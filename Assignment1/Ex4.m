clc
clear
close all;
%% code parameters

n = 7; k = 4; r = n - k;

I = eye(k);
P = [   1 1 0;
        1 0 1;
        0 1 1;
        1 1 1   ];
G = [I P];


%% verify that the code is perfect

% generate all possible information vectors
V = dec2bin(0:2^k-1) - '0';

% compute the codebook
C = mod(V * G, 2);

% compute Hamming weights
w = sum(C, 2);

% compute weight multiplicities A_i for i = 0, ..., n
% the i-th element corresponds to weight i => A is indexed with the weight
A = zeros(1, n + 1); % preallocation

for i = 0:n
    A(i+1) = sum(w == i); % sum of all the w1 elements that have value = i
end

% minimum distances and their multiplicities
nonzero_w = w(w > 0); % only consider non-zero weights

% find the smallest non-zero weight
dmin = min(nonzero_w); 

% find the multiplicity corresponding to the smallest weight
Amin = A(dmin+1);

t = floor((dmin-1)/2);

% compute all possible received vectors
Y = dec2bin(0:2^n-1) - '0';


% compute decision regions
D_c = cell(size(C, 1), 1);

% for each received vector, calculate Hamming distance to all codewords
for j = 1:size(Y, 1)
    y = Y(j, :);

    % compute Hamming distances from y to all codewords
    d_H = zeros(size(C, 1), 1);
    for i = 1:size(C, 1)
        d_H(i) = sum(y ~= C(i, :));
    end

    % find the minimum Hamming distance
    d_min = min(d_H);

    % find all codewords achieving the minimum distance
    d_min_idx = find(d_H == d_min);

    % include received vector in the decision region of each matching codeword;
    % in case of a tie, the received vector is included in all corresponding decision regions
    for idx = d_min_idx'
        D_c{idx} = [D_c{idx}; y];
    end
end


% compute Hamming spheres
S_c = cell(size(C, 1), 1);

% for each received vector, calculate Hamming distance to all codewords
for j = 1:size(Y, 1)
    y = Y(j, :);

    % compute Hamming distances from y to all codewords
    d_H = zeros(size(C, 1), 1);
    for i = 1:size(C, 1)
        d_H(i) = sum(y ~= C(i, :));
    end

    % find all codewords less or equal to t
    d_idx = find(d_H <= t);

    % include received vector in the decision region of each matching codeword;
    % in case of a tie, the received vector is included in all corresponding decision regions
    for idx = d_idx'
        S_c{idx} = [S_c{idx}; y];
    end
end

perfect = true;

% a code is perfect is the Hamming spheres of radius t exactly cover the
% received-vector space without overlap
for i = 1:size(C, 1)
    perfect = perfect && isequal(D_c{i}, S_c{i});
end

if perfect
    fprintf('VERIFICATION: The (7,4) Hamming code is PERFECT.\n\n');
else
    fprintf('VERIFICATION: The (7,4) code is NOT perfect.\n\n');
end

%% word error probability
p = logspace(-6, -0.5, 1000);
P_c = zeros(size(p)); % probability of correct decoding

for i = 0:t
    P_c = P_c + nchoosek(n,i) * p.^i .* (1-p).^(n-i);
end

P_w = 1 - P_c; % probability of a codeword error
P_w_perfect = P_w;

figure;
axis tight
loglog(p, P_w, 'LineWidth', 1.5); grid on;
xlabel('$p$', 'Interpreter','latex')
ylabel('$P_w(E)$', 'Interpreter','latex')
% title('Word error probability');
set(gca, "FontSize", 20)

%% Monte Carlo simulation
p_MC   = logspace(-2, -1, 10);
target = 100;
alpha  = 0.05;


% Initialization
P_W_MC  = zeros(1, length(p_MC));
CI_low  = zeros(1, length(p_MC));
CI_high = zeros(1, length(p_MC));


rng(12); % reproducibility

for i = 1:length(p_MC)
    pc = p_MC(i);
    N_WE = 0; N_trans = 0;
    while N_WE < target
        v = randi([0 1], 1, k);
        c = mod(v * G, 2);
        e = rand(1, n) < pc;
        y = xor(c, e);

        d_H = sum(C ~= y, 2);
        [~, decoded_idx] = min(d_H);

        v_hat = V(decoded_idx, :);

        N_WE    = N_WE + any(v_hat ~= v);
        N_trans = N_trans + 1;
    end
    P_W_MC(i) = N_WE / N_trans;

    r = N_WE; f = N_trans - N_WE;
    CI_low(i)  = betainv(alpha/2,   r,   f+1);
    CI_high(i) = betainv(1-alpha/2, r+1, f);
end

% analytical curve
p_th = logspace(-2, -1, 200);

P_c = zeros(size(p_th)); % probability of correct decoding

for i = 0:t
    P_c = P_c + nchoosek(n,i) * p_th.^i .* (1-p_th).^(n-i);
end

P_th = 1 - P_c; 

% plot with confidence intervals
figure;
loglog(p_th, P_th, 'b-', 'LineWidth', 1.5); hold on;
%loglog(p_MC, P_W_theory_MC, 'ks', 'MarkerFaceColor', 'k'); % theory at MC points
for j = 1:length(target)
    errorbar(p_MC, P_W_MC, P_W_MC-CI_low, CI_high-P_W_MC, '-o', 'LineWidth', 1.2);
end
set(gca, 'XScale', 'log', 'YScale', 'log', 'FontSize', 15);
grid on;
xticks(p_MC);
axis tight
xlabel('BSC crossover probability p'); ylabel('P_W(E)');
legend('Analytical', 'MC: 100 word errors', 'Location', 'northwest');


%% Another (8,4) code

n = 8; k = 4; r = n - k;

I = eye(4);
P = [   1 1 0 1;
        1 0 1 1;
        0 1 1 1;
        1 1 1 0  ];
G = [I P];

%% verify that the code is perfect

% generate all possible information vectors
V = dec2bin(0:2^k-1) - '0';

% compute the codebook
C = mod(V * G, 2);

% compute Hamming weights
w = sum(C, 2);

% compute weight multiplicities A_i for i = 0, ..., n
% the i-th element corresponds to weight i => A is indexed with the weight
A = zeros(1, n + 1); % preallocation

for i = 0:n
    A(i+1) = sum(w == i); % sum of all the w1 elements that have value = i
end

% minimum distances and their multiplicities
nonzero_w = w(w > 0); % only consider non-zero weights

% find the smallest non-zero weight
dmin = min(nonzero_w); 

% find the multiplicity corresponding to the smallest weight
Amin = A(dmin+1);

t = floor((dmin-1)/2);

% compute all possible received vectors
Y = dec2bin(0:2^n-1) - '0';


% compute decision regions
D_c = cell(size(C, 1), 1);

% for each received vector, calculate Hamming distance to all codewords
for j = 1:size(Y, 1)
    y = Y(j, :);

    % compute Hamming distances from y to all codewords
    d_H = zeros(size(C, 1), 1);
    for i = 1:size(C, 1)
        d_H(i) = sum(y ~= C(i, :));
    end

    % find the minimum Hamming distance
    d_min = min(d_H);

    % find all codewords achieving the minimum distance
    d_min_idx = find(d_H == d_min);

    % include received vector in the decision region of each matching codeword;
    % in case of a tie, the received vector is included in all corresponding decision regions
    for idx = d_min_idx'
        D_c{idx} = [D_c{idx}; y];
    end
end


% compute Hamming spheres
S_c = cell(size(C, 1), 1);

% for each received vector, calculate Hamming distance to all codewords
for j = 1:size(Y, 1)
    y = Y(j, :);

    % compute Hamming distances from y to all codewords
    d_H = zeros(size(C, 1), 1);
    for i = 1:size(C, 1)
        d_H(i) = sum(y ~= C(i, :));
    end

    % find all codewords less or equal to t
    d_idx = find(d_H <= t);

    % include received vector in the decision region of each matching codeword;
    % in case of a tie, the received vector is included in all corresponding decision regions
    for idx = d_idx'
        S_c{idx} = [S_c{idx}; y];
    end
end

perfect = true;

% a code is perfect is the Hamming spheres of radius t exactly cover the
% received-vector space without overlap
for i = 1:size(C, 1)
    perfect = perfect && isequal(D_c{i}, S_c{i});
end

if perfect
    fprintf('VERIFICATION: The (8,4) Hamming code is PERFECT.\n\n');
else
    fprintf('VERIFICATION: The (8,4) code is NOT perfect.\n\n');
end

%% Analytical error codeword probability
E   = dec2bin(0:2^n-1, n) - '0'; % all error patterns (y = e, all-zero word sent)
w_e = sum(E, 2); % Hamming weight = distance from the zero codeword

% D(0): vectors y whose distance from 0 is minimal among all codewords
in_D0 = false(2^n, 1);
for m = 1:2^n
    d = sum(C ~= E(m,:), 2); % distances from y to all codewords
    in_D0(m) = (w_e(m) == min(d)); % C(1,:) is the zero codeword, so d(1) = w_e(m)
end
D0_idx = find(in_D0);

p   = 0:0.01:1;
P_c = zeros(size(p)); % probability of correct decoding
for i = 1:numel(D0_idx)
    m   = D0_idx(i);
    P_c = P_c + p.^w_e(m) .* (1-p).^(n-w_e(m));
end
P_w = 1 - P_c; % probability of a word error

P_w_not_perfect = P_w;

figure;
plot(p, P_w, 'LineWidth', 1.5); grid on;
xlabel('$p$', 'Interpreter','latex')
ylabel('$P_w(p)$', 'Interpreter','latex')
title('Word error probability');
set(gca, 'FontSize', 20)

%% Monte Carlo simulation

p_MC   = logspace(-2, -1, 10);
target = 100;
alpha  = 0.05;


% Initialization
P_W_MC  = zeros(1, length(p_MC));
CI_low  = zeros(1, length(p_MC));
CI_high = zeros(1, length(p_MC));


rng(12); % reproducibility

for i = 1:length(p_MC)
    pc = p_MC(i);
    N_WE = 0; N_trans = 0;
    while N_WE < target
        v = randi([0 1], 1, k);
        c = mod(v * G, 2);
        e = rand(1, n) < pc;
        y = xor(c, e);

        d_H = sum(C ~= y, 2);
        idx_tx = v * (2.^(k-1:-1:0))' + 1; % row of C of the transmitted word

        % error only if some other codeword is strictly closer
        N_WE    = N_WE + (d_H(idx_tx) > min(d_H));
        N_trans = N_trans + 1;
    end
    P_W_MC(i) = N_WE / N_trans;

    r = N_WE; f = N_trans - N_WE;
    CI_low(i)  = betainv(alpha/2, r,f+1);
    CI_high(i) = betainv(1-alpha/2, r+1, f);
end

% analytical curve
p_th = logspace(-2, -1, 200);
P_c = zeros(size(p_th));                   % probability of correct decoding
for i = 1:numel(D0_idx)
    m   = D0_idx(i);
    P_c = P_c + p_th.^w_e(m) .* (1-p_th).^(n-w_e(m));
end
P_th = 1 - P_c;                          % probability of a word error

% plot with confidence intervals
figure;
loglog(p_th, P_th, 'b-', 'LineWidth', 1.5); hold on;
for j = 1:length(target)
    errorbar(p_MC, P_W_MC, P_W_MC-CI_low, CI_high-P_W_MC, '-o', 'LineWidth', 1.2);
end
set(gca, 'XScale', 'log', 'YScale', 'log', 'FontSize', 15);
grid on;
xticks(p_MC);
axis tight
xlabel('BSC crossover probability p'); ylabel('P_w(E)');
legend('Analytical', 'MC: 100 word errors', 'Location', 'northwest');


%% Compare the two codes
p = 0:0.01:1;

% Re-evaluate C(7,4) over p = 0:0.01:1
P_c_74 = (1-p).^7 + 7.*p.*(1-p).^6;
P_w_perfect = 1 - P_c_74;

% Re-evaluate C(8,4) over p = 0:0.01:1
P_c_84 = zeros(size(p));
for i = 1:numel(D0_idx)
    m = D0_idx(i);
    P_c_84 = P_c_84 + p.^w_e(m) .* (1-p).^(n-w_e(m));
end
P_w_not_perfect = 1 - P_c_84;

figure;
plot(p, P_w_perfect, p, P_w_not_perfect, 'LineWidth', 1.5);
grid on;
xlabel('$p$', 'Interpreter','latex');
ylabel('$P_w(E)$', 'Interpreter','latex');
title('Word Error Probability Comparison');
legend('$P_w(E)$ for $\mathcal{C}(7,4)$', '$P_w(E)$ for $\mathcal{C}(8,4)$', 'Interpreter', 'latex', 'Location', 'southeast');
set(gca, 'FontSize', 18);