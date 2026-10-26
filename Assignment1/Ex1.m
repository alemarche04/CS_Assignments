clear; clc; close all;
rng(2); % reproducibility

n = 10;  k = 5;

% generator matrices from the assignment
I5 = eye(k);
P1 = [1 1 1 0 0;
      0 1 1 1 0;
      0 0 1 1 1;
      1 0 0 1 1;
      1 1 0 0 1];
P2 = [1 1 0 0 0;
      0 1 1 1 0;
      0 0 1 1 1;
      1 0 0 1 1;
      1 1 0 0 1];
G1 = [I5 P1];
G2 = [I5 P2];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part1: Codebooks and Weight Spectra 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% generate all 2^k = 32 information vectors
% dec2bin - Convert decimal integer to its binary representation
V = dec2bin(0:2^k-1) - '0'; % - '0' => conversion from binary strings to binary values

% encode all information vectors
C1 = mod(V * G1, 2);
fprintf('\nG1 codebook:\n');
disp(C1);

C2 = mod(V * G2, 2);
fprintf('\nG2 codebook:\n');
disp(C2);

% Hamming weights of all codewords
w1 = sum(C1, 2);
w2 = sum(C2, 2);

% compute weight multiplicities A_i for i = 0..n
% the i-th element corresponds to weight i => A is indexed with the weight
A1 = zeros(1, n + 1); % preallocation
A2 = zeros(1, n + 1);

for i = 0:n
    A1(i+1) = sum(w1 == i); % sum of all the w1 elements that have value = i
    A2(i+1) = sum(w2 == i);
end

fprintf('G1 multiplicities:\n');
for j = 1:length(A1)
    fprintf('A%i = %i\n', j-1, A1(j))
end

fprintf('G2 multiplicities:\n');
for j = 1:length(A2)
    fprintf('A%i = %i\n', j-1, A2(j))
end


% minimum distances and their multiplicities
nonzero_w1 = w1(w1 > 0); % only consider non-zero weights
fprintf('G1 weights:\n');
disp(w1.');

nonzero_w2 = w2(w2 > 0);
fprintf('G2 weights:\n');
disp(w2.');

dmin1 = min(nonzero_w1); % find the smallest non-zero weight
dmin2 = min(nonzero_w2);

Amin1 = A1(dmin1+1); % find the multiplicity corresponding to the smallest weight
Amin2 = A2(dmin2+1);


fprintf('G1: d_min = %d, A_min = %d\n', dmin1, Amin1);

fprintf('G2: d_min = %d, A_min = %d\n', dmin2, Amin2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part 2: Analytical Performance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% P_UE(p) = sum_{j=1}^n A_j p^j (1-p)^(n-j)
% NOTE: A_0 is not included because the all-zero error corresponds to correct reception

p_ax = logspace(-6, -1, 1000); % p_ax ​= [10^−6, ...,10^−1] with 1000 points

P_UE1 = zeros(size(p_ax)); % preallocation
P_UE2 = zeros(size(p_ax));

for m = 1:length(p_ax) % implement P_UE(p) = sum_{j=1}^n A_j p^j (1-p)^(n-j)
    p = p_ax(m);
    for j = 1:n
        P_UE1(m) = P_UE1(m) + A1(j+1) * p^j * (1-p)^(n-j);
        P_UE2(m) = P_UE2(m) + A2(j+1) * p^j * (1-p)^(n-j);
    end
end

% plot
figure;
loglog(p_ax, P_UE1, 'LineWidth', 1.5);
hold on;
loglog(p_ax, P_UE2, 'LineWidth', 1.5);

grid on;
xlabel('BSC crossover probability p'); ylabel('P_{UE}(p)');
legend('Code G_1','Code G_2','Location','northwest');
set(gca, 'FontSize', 15);


% Leading terms for small p:
% Code 1: A_min * p^dmin
% Code 2: A_min * p^dmin

% leading-term approximations for small p
PUE1_approx = Amin1 * p_ax.^dmin1;
PUE2_approx = Amin2 * p_ax.^dmin2;

% plot analytical curves and leading-term approximations
figure;
loglog(p_ax, P_UE1, 'b-', 'LineWidth', 1.5);
hold on;
loglog(p_ax, PUE1_approx, 'c--', 'LineWidth', 1.5);
loglog(p_ax, P_UE2, 'r-', 'LineWidth', 1.5);
loglog(p_ax, PUE2_approx, 'm--', 'LineWidth', 1.5);

grid on;
xlabel('BSC crossover probability p');
ylabel('P_{UE}(p)');
set(gca, 'FontSize', 15);
legend('Code G_1', 'G_1: A_{1, min} * p^{d_{1, min}}', 'Code G_2', 'G_2: A_{2, min} * p^{d_{2, min}}', 'Location', 'northwest');

fprintf('\nG1 leading term for small p: PUE1 ~= %d * p^%d\n', Amin1,dmin1);
fprintf('G2 leading term for small p: PUE2 ~= %d * p^%d\n', Amin2,dmin2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part 3: Syndrome detector (code G1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% For a systematic generator matrix G = [I_k P] the parity check matrix is  
% H^T = [      P 
%           I_(n-k)     ]
% The syndrome of a row vector y is s = y H^T
r = n-k;

HT = [P1; eye(r)];
fprintf('\nH^T = \n'); disp(HT)

% NOTE: y*H^T = 0 IIF y is a codeword

% choose received vector y1 that is a codeword
y1 = C1(2,:);
fprintf('\ny1 = '); disp(y1)

% y1 syndrome (y1 is a codeword => s1 = 0)
s1 = mod(y1 * HT, 2); 
fprintf('\nsyndrome s1 = y1 * H^T = '); disp(s1)

% choose a received vector y2 = codeword + one-bit error
y2 = C1(2,:);
error_position = 3;
y2(error_position) = mod(y2(error_position) + 1, 2);
fprintf('\ny2 = '); disp(y2)

% y2 syndrome (y2 is not a codeword => s2 = 0)
s2 = mod(y2 * HT, 2); 
fprintf('\nsyndrome s2 = y2 * H^T = '); disp(s2)


% obtain a different parity check matrix by swapping two or more columns of
% the original matrix HT

% swap the first two columns
HT1 = HT;
HT1(:, 1)    = HT(:, 2);
HT1(:, 2)    = HT(:, 1);

% a different (10x5) parity-check matrix can be obtained by multiplying H^T by an invertible (5x5) matrix
% (any invertible transformation preserves the condition s = 0)

M =  [1 0 0 0 0; % invertible matrix
      0 0 0 1 0;
      0 0 1 0 0;
      0 0 0 0 1;
      0 1 0 0 0];

H1T = mod(HT*M, 2); % matrix multiplication between H^T and M
fprintf('\nH1^T = \n'); disp(H1T)

s1_H1 = mod(y1 * H1T, 2); % y1 is a codeword => s1_H1 = 0
s2_H1 = mod(y2 * H1T, 2); % y2 is a not codeword => s2_H1 != 0

% y1 syndrome
fprintf('\ny1 = '); disp(y1)
fprintf('\nsyndrome s1 with new matrix H1^T = '); disp(s1_H1)

% y2 syndrome
fprintf('\ny2 = '); disp(y2)
fprintf('\nsyndrome s2 with new matrix H1^T = '); disp(s2_H1)


% a new matrix H2^T with dimensions (10x7) is still a parity-check matrix if we append two linearly dependent columns
% for example: (col6 = col1 + col2,  col7 = col4 + col5) => rank stays 5

H2T = [HT (mod(HT(:, 1) + HT(:, 2), 2)) (mod(HT(:, 4) + HT(:, 5), 2))];
fprintf('\nH2^T = \n'); disp(H2T)

s1_H2 = mod(y1 * H2T, 2); % y1 is a codeword => s1_H2 = 0
s2_H2 = mod(y2 * H2T, 2); % y2 is a not codeword => s2_H2 != 0

% y1 syndrome
fprintf('\ny1 = '); disp(y1)
fprintf('\nsyndrome s1 with new matrix H2^T = '); disp(s1_H2)

% y2 syndrome
fprintf('\ny2 = '); disp(y2)
fprintf('\nsyndrome s2 with new matrix H2^T = '); disp(s2_H2)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part 4: Monte Carlo Simulation (code G1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p_MC = [0.10, 0.09, 0.08, 0.07, 0.06];

target = [10 100]; % stopping criteria

% initialization
P_UE_MC = zeros(length(target), length(p_MC));

CI_low  = zeros(length(target), length(p_MC));
CI_high = zeros(length(target), length(p_MC));

for j = 1:length(target)
    for i = 1:length(p_MC)
        p = p_MC(i);

        N_UE = 0;
        N_trans = 0;


        while N_UE < target(j)

            v = randi([0 1], 1, k); % random information vector of length k
            c = mod(v * G1, 2); % compute codeword from random information vector

            % transmission over BSC
            e = rand(1, n) < p; % error vector: if i-th value of rand is < p => e(i) = 1

            y = mod(c + e, 2); % received vector

            s = mod(y*HT, 2); % syndrome of received vector

            % there is an undetected error if we have zero syndrome and y differs from c
            if all(s == 0) && any(y ~= c)
                N_UE = N_UE + 1;
            end

            N_trans = N_trans + 1;
        end

        P_UE_MC(j, i) = N_UE / N_trans;

        % 95% confidence interval for the stopping rule
        alpha = 0.05;

        r = N_UE;
        f = N_trans - N_UE;

        % exact Clopper-Pearson 95% confidence bounds (alpha = 0.05) are 
        % calculated using the inverse Beta cumulative distribution 
        % function (betainv)
        CI_low(j, i)  = betainv(alpha/2, r, f + 1);
        CI_high(j, i) = betainv(1-alpha/2, r, f);

    end
end

% analytical curve from Part 2 => P_UE1

% plot comparison
figure;
loglog(p_ax, P_UE1, 'b-', 'LineWidth', 1.5);
hold on;
for j = 1:length(target)
    loglog(p_MC, P_UE_MC(j,:), 'o-', 'LineWidth', 1.5);
end
grid on;
xlabel('BSC crossover probability p'); ylabel('P_{UE}(p)');
legend('Analytical', 'Monte Carlo: 10 UE', 'Monte Carlo: 100 UE', 'Location', 'northwest');
xlim([min(p_MC), max(p_MC)]);
set(gca, 'FontSize', 15);


figure;
loglog(p_ax, P_UE1, 'b-', 'LineWidth', 1.5);
hold on;

for j = 1:length(target)

    lower_error = P_UE_MC(j,:) - CI_low(j,:);
    upper_error = CI_high(j,:) - P_UE_MC(j,:);

    errorbar(p_MC, P_UE_MC(j,:), ...
        lower_error, upper_error, 'o-', ...
        'LineWidth', 1.2);

end

grid on;
xlabel('BSC crossover probability p');
ylabel('P_{UE}(p)');

legend('Analytical', ...
       'Monte Carlo: 10 UE', ...
       'Monte Carlo: 100 UE', ...
       'Location', 'northwest');

xlim([min(p_MC), max(p_MC)]);
set(gca, 'FontSize', 15);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part 5: Retransmissions (ARQ) (Code G1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% For one transmission:
%
% p0 = P(COR, i=1)
%    = probability of no channel errors
%    = (1-p)^n
%
% p1 = P(UE, i=1)
%    = probability of undetected error
%    = sum_{i=d_min}^{n} A_i * (1-p)^(n-i) * p^i
%    = P_UE(p)
%
% p2 = P(REJ, i=1)
%    = probability of rejection
%    = 1 - p0 - p1

% With at most i attempts, correct reception occurs if the receiver
% gets an error-free codeword (syndrome = 0) after 0,1,...,i-1 rejected attempts
% P(COR,i) = p0 * sum_{j=0}^{i-1} p2^j

p_ax_ARQ = logspace(-2, -1, 1000); % p_ax_ARQ ​= [10^−2, ...,10^−1] with 1000 points

% initialization
p0_ARQ = zeros(size(p_ax_ARQ));
p1_ARQ = zeros(size(p_ax_ARQ));
p2_ARQ = zeros(size(p_ax_ARQ));

for m = 1:length(p_ax_ARQ) % compute p0, p1 and p2 for each p_ax_ARQ ​= [10^−2, ...,10^−1]
    p = p_ax_ARQ(m);

    p0_ARQ(m) = (1-p)^n;

    p1_ARQ(m) = 0;
    for j = 1:n % implement the formula sum_{i=d_min}^{n} A_i * (1-p)^(n-i) * p^i
        p1_ARQ(m) = p1_ARQ(m) + ...
            A1(j+1)*p^j*(1-p)^(n-j);
    end

    p2_ARQ(m) = 1 - p0_ARQ(m) - p1_ARQ(m);
end

P_COR_1 = p0_ARQ;
P_COR_2 = p0_ARQ .* (1 + p2_ARQ);
P_COR_3 = p0_ARQ .* (1 + p2_ARQ + p2_ARQ.^2);

figure;
loglog(p_ax_ARQ, 1 - P_COR_1, 'LineWidth', 1.5);
hold on;
semilogx(p_ax_ARQ, 1 - P_COR_2, 'LineWidth', 1.5);
semilogx(p_ax_ARQ, 1 - P_COR_3, 'LineWidth', 1.5);

grid on;
xlabel('BSC crossover probability p'); ylabel('1 - P(COR, i) (error probability)');
legend('i = 1','i = 2','i = 3','Location','northwest');
set(gca, 'FontSize', 15);