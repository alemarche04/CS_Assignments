clc
clear
close all;
%% code parameters
m = 6;
n = 2^m -1; k = 57; r = n - k;

t = 1; % error correction capability

% Systematic Hamming (63, 57) matrices
% H is (6 x 63), G is (57 x 63)
[H, G] = hammgen(m);

HT = H';

%% word error probability
p = logspace(-5, -1, 100);
P_c = zeros(size(p)); % probability of correct decoding

for i = 0:t
    P_c = P_c + nchoosek(n,i) * p.^i .* (1-p).^(n-i);
end

P_w = 1 - P_c; % probability of a codeword error

% P_w = 1 - (1 - p).^63 - 63*p.*(1 - p).^62;

figure;
loglog(p, P_w, 'LineWidth', 1.5); grid on;
xlabel('$p$', 'Interpreter','latex')
ylabel('$P_w(E)$', 'Interpreter','latex')
% title('Word error probability');
axis tight
set(gca, "FontSize", 20)


%% bounded-distance syndrome decoder

% construct Syndrome Lookup Table for Radius t = 1
% map 6-bit syndrome vector with corresponding 63-bit error vector e
syndrome_table = containers.Map('KeyType', 'char', 'ValueType', 'any');

% all-zero syndrome => all-zero error vector
syndrome_table(char('0' + zeros(1, r))) = zeros(1, n);

% Single-bit error patterns (63 entries)
% errors that can be surely corrected (t=1)
for i = 1:n
    e_single = zeros(1, n);
    e_single(i) = 1;
    s = mod(e_single * HT, 2); % syndrome s = e * H^T
    syndrome_table(char('0' + s)) = e_single;
end

% analytical Codeword Error Probability
p_vec = logspace(-2, -0.5, 10); % range of channel bit error probabilities p
P_w_analytical = zeros(size(p_vec));

for idx = 1:length(p_vec)
    p = p_vec(idx);
    P_correct = 0;
    for i = 0:t
        P_correct = P_correct + nchoosek(n, i) * (p^i) * ((1 - p)^(n - i));
    end
    P_w_analytical(idx) = 1 - P_correct;
end

% Monte Carlo simulation of Bounded-Distance Syndrome Decoder
rng(22); % reproducibility
target_errors = 100;
alpha = 0.05; % parameter for confidence intervals

P_w_sim     = zeros(size(p_vec));
CI_low1     = zeros(1, length(p_vec));
CI_high1    = zeros(1, length(p_vec));

for j = 1:length(p_vec)
    p = p_vec(j);
    word_errors = 0;
    N_WE = 0; N_trans = 0;
    while N_WE < target_errors
        % transmit zero codeword
        c_tx = zeros(1, n); 

        % BSC Channel Noise: bit flips with probability p
        e_channel = rand(1, n) < p;
        y = mod(c_tx + e_channel, 2); % Received vector y = c + e

        % Bounded-Distance Syndrome Decoder
        s = mod(y * HT, 2); % compute syndrome s = y * H^T
        s_key = char('0' + s);

        % look up estimated error vector e_hat
        if isKey(syndrome_table, s_key)
            e_hat = syndrome_table(s_key);
        else
            % syndrome not in table (occurs when > t errors happen)
            e_hat = zeros(1, n); % decoding failure
        end

        % correct received vector c_hat = y + e_hat
        c_hat = mod(y + e_hat, 2);

        % check for block decoding error
        if any(c_hat ~= c_tx)
            N_WE = N_WE + 1;
        end
        N_trans = N_trans + 1;
        
    end

    P_w_sim(j) = N_WE / N_trans;

    r = N_WE; f = N_trans - N_WE;
    CI_low1(j)  = betainv(alpha/2, r, f+1);
    CI_high1(j) = betainv(1-alpha/2, r+1, f);
end

% plotting & Performance Comparison
figure;
loglog(p_vec, P_w_analytical, 'b-', 'LineWidth', 2);
hold on
errorbar(p_vec, P_w_sim, P_w_sim - CI_low1, CI_high1 - P_w_sim, '-ro', ...
    'LineWidth', 1.5, 'MarkerSize', 8, 'MarkerFaceColor', 'r');
hold off
grid on; axis tight;
xlabel('BSC crossover probability p', 'Interpreter', 'latex');
ylabel('$P_w(E)$', 'Interpreter', 'latex');
%title('Hamming (63,57) Syndrome Decoding Performance');
legend('Analytical $P_w(E)$', 'Monte Carlo (100 word errors)', 'Interpreter', 'latex', 'Location', 'southeast');
set(gca, 'FontSize', 18);