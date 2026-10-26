%% ========================================================================
%  Communication Systems - Assignment 1 - Exercise 1
%  Error detection with two (10,5) systematic linear block codes
%  Requires: base MATLAB only (no toolboxes). R2016b+ (local functions in scripts).
%  ========================================================================
clear; clc; close all;
rng(2026);                         % reproducibility

n = 10;  k = 5;

%% ------------------------------------------------------------------------
%  Generator matrices  G = [I5 | P]
%  ------------------------------------------------------------------------
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

%% ========================================================================
%  PART 1 - Codebooks and weight spectra
%  ========================================================================
[C1, w1, A1, dmin1, Amin1] = codebook(G1);
[C2, w2, A2, dmin2, Amin2] = codebook(G2);

fprintf('=== PART 1: weight spectra ===\n');
fprintf('i        :'); fprintf('%4d', 0:n); fprintf('\n');
fprintf('A_i  (C1):'); fprintf('%4d', A1);  fprintf('   | sum = %d\n', sum(A1));
fprintf('A_i  (C2):'); fprintf('%4d', A2);  fprintf('   | sum = %d\n', sum(A2));
fprintf('Code 1: dmin = %d, Amin = %d\n', dmin1, Amin1);
fprintf('Code 2: dmin = %d, Amin = %d\n\n', dmin2, Amin2);

figure('Name','Weight spectra');
bar(0:n, [A1(:) A2(:)]); grid on;
xlabel('Hamming weight i'); ylabel('A_i');
legend('Code 1 (G_1)','Code 2 (G_2)','Location','northwest');
title('Weight spectra of the two codes');

%% ========================================================================
%  PART 2 - Analytical P_UE(p)
%  P_UE(p) = sum_{i=1}^{n} A_i p^i (1-p)^(n-i)
%  ========================================================================
p = logspace(-6, -1, 200);
PUE1 = PUE(A1, p);
PUE2 = PUE(A2, p);

% Leading-term approximations for p -> 0:  P_UE ~ Amin * p^dmin
lead1 = Amin1 * p.^dmin1;
lead2 = Amin2 * p.^dmin2;

figure('Name','Analytical PUE');
loglog(p, PUE1, 'b-',  'LineWidth', 1.8); hold on;
loglog(p, PUE2, 'r-',  'LineWidth', 1.8);
loglog(p, lead1, 'b--', 'LineWidth', 1.0);
loglog(p, lead2, 'r--', 'LineWidth', 1.0);
grid on; xlabel('p'); ylabel('P_{UE}(p)');
legend('Code 1','Code 2', ...
       sprintf('A_{min} p^{d_{min}} (C1: %d p^{%d})', Amin1, dmin1), ...
       sprintf('A_{min} p^{d_{min}} (C2: %d p^{%d})', Amin2, dmin2), ...
       'Location','southeast');
title('Probability of undetected error over the BSC');

fprintf('=== PART 2: leading terms ===\n');
fprintf('Code 1: PUE ~ %d * p^%d\n', Amin1, dmin1);
fprintf('Code 2: PUE ~ %d * p^%d\n\n', Amin2, dmin2);

%% ========================================================================
%  PART 3 - Syndrome detector (code 1)
%  ========================================================================
fprintf('=== PART 3: syndrome detector ===\n');

% 3.1  Parity-check matrix, lecture convention: H^T = [P; I5]  (10x5)
HT = [P1; I5];
assert(all(all(mod(G1*HT,2)==0)), 'G*H^T must be zero');

% 3.2  Two received vectors
v   = [1 0 1 1 0];
c   = mod(v*G1, 2);                 % a codeword  -> zero syndrome
y1  = c;
y2  = c;  y2(3) = 1 - y2(3);        % single bit flip -> non-zero syndrome
s1  = mod(y1*HT, 2);
s2  = mod(y2*HT, 2);
fprintf('y1 = %s  (codeword)           -> s = %s\n', mat2str(y1), mat2str(s1));
fprintf('y2 = %s  (1 bit flipped)      -> s = %s\n', mat2str(y2), mat2str(s2));

% 3.3  Alternative H1^T (10x5): H1^T = H^T * M, with M invertible over F2.
%      Same column space => same null space (left kernel) => same code.
M  = [1 1 0 0 0;
      0 1 1 0 0;
      0 0 1 1 0;
      0 0 0 1 1;
      0 0 0 0 1];                   % upper bidiagonal, det = 1 (invertible)
H1T = mod(HT*M, 2);
assert(verifyH(H1T, C1, n), 'H1T is not a valid parity-check matrix');
fprintf('H1^T (10x5) valid. s(y1) = %s, s(y2) = %s\n', ...
        mat2str(mod(y1*H1T,2)), mat2str(mod(y2*H1T,2)));

% 3.4  Redundant H2^T (10x7): append two linearly dependent columns
%      (col6 = col1 + col2,  col7 = col4 + col5), rank stays 5.
Ext = [eye(5), [1;1;0;0;0], [0;0;0;1;1]];
H2T = mod(HT*Ext, 2);
assert(verifyH(H2T, C1, n), 'H2T is not a valid parity-check matrix');
fprintf('H2^T (10x7) valid. s(y1) = %s, s(y2) = %s\n\n', ...
        mat2str(mod(y1*H2T,2)), mat2str(mod(y2*H2T,2)));

%% ========================================================================
%  PART 4 - Monte Carlo simulation (code 1)
%  ========================================================================
fprintf('=== PART 4: Monte Carlo ===\n');
pMC     = [0.10 0.09 0.08 0.07 0.06];
targets = [10 100];                 % stopping criteria (undetected errors)
batch   = 2e5;                      % vectorised block size
z       = 1.96;                     % 95% confidence level

PUE_MC = zeros(numel(targets), numel(pMC));
Ntx    = zeros(numel(targets), numel(pMC));
CIlo   = zeros(size(PUE_MC));
CIhi   = zeros(size(PUE_MC));

for it = 1:numel(targets)
    for ip = 1:numel(pMC)
        tic;
        [nUE, N] = mc_until(G1, HT, pMC(ip), targets(it), batch);
        PUE_MC(it,ip) = nUE / N;
        Ntx(it,ip)    = N;
        % Normal-approximation (Wald) 95% confidence interval
        half = z * sqrt(PUE_MC(it,ip)*(1-PUE_MC(it,ip)) / N);
        CIlo(it,ip) = max(PUE_MC(it,ip) - half, eps);
        CIhi(it,ip) = PUE_MC(it,ip) + half;
        fprintf('target=%3d  p=%.2f  N=%9d  PUE_MC=%.3e  PUE_th=%.3e  (%.1fs)\n', ...
                targets(it), pMC(ip), N, PUE_MC(it,ip), PUE(A1,pMC(ip)), toc);
    end
end

pFine = linspace(0.06, 0.10, 100);
figure('Name','Monte Carlo vs analytical');
semilogy(pFine, PUE(A1,pFine), 'k-', 'LineWidth', 1.8); hold on;
cols = {'b','r'}; mk = {'o','s'};
for it = 1:numel(targets)
    errorbar(pMC, PUE_MC(it,:), PUE_MC(it,:)-CIlo(it,:), CIhi(it,:)-PUE_MC(it,:), ...
             [cols{it} mk{it} '--'], 'LineWidth', 1.2, 'MarkerSize', 7);
end
grid on; xlabel('p'); ylabel('P_{UE}');
legend('Analytical','MC (10 undetected errors)','MC (100 undetected errors)', ...
       'Location','southeast');
title('Monte Carlo estimates with 95% confidence intervals (code 1)');

figure('Name','Simulation cost');
semilogy(pMC, Ntx(1,:), 'bo-', pMC, Ntx(2,:), 'rs-', 'LineWidth', 1.2); grid on;
xlabel('p'); ylabel('Transmitted codewords N');
legend('10 errors','100 errors'); title('Computational cost of the simulation');

%% ========================================================================
%  PART 5 - Ideal error-detection ARQ (code 1)
%  ========================================================================
%  Per-attempt probabilities:
%   p0 = P(COR,1) = (1-p)^n                       (error pattern = 0)
%   p1 = P(Ue,1)  = sum_{i>=1} A_i p^i (1-p)^(n-i) = PUE(p)
%   p2 = P(REJ,1) = 1 - sum_{i>=0} A_i p^i (1-p)^(n-i)
%  With at most i attempts:
%   P(COR,i) = p0 * sum_{j=0}^{i-1} p2^j = p0 (1 - p2^i)/(1 - p2)
%   P(Ue,i)  = p1 * sum_{j=0}^{i-1} p2^j
%   P(fail,i)= p2^i
pA = logspace(-2, -1, 100);
p0 = (1-pA).^n;
p1 = PUE(A1, pA);
p2 = 1 - (p0 + p1);

figure('Name','ARQ');
styles = {'b-','r-','m-'};
for i = 1:3
    Pcor = p0 .* geomsum(p2, i);
    loglog(pA, 1 - Pcor, styles{i}, 'LineWidth', 1.8); hold on;
end
grid on; xlabel('p'); ylabel('1 - P(COR, i)');
legend('i = 1','i = 2','i = 3','Location','southeast');
title('Probability of not receiving the correct codeword (ARQ)');

%% ========================================================================
%  PART 6 - Expected number of transmissions and conditional P(UE)
%  ========================================================================
%  T = number of transmissions with max i attempts:
%   P(T=j) = p2^(j-1) (1-p2)  for j<i ;   P(T=i) = p2^(i-1)
%   E[T_i] = sum_{j=0}^{i-1} p2^j = (1 - p2^i)/(1 - p2)  -> 1/(1-p2) as i->inf
figure('Name','Expected transmissions');
for i = 1:3
    ET = geomsum(p2, i);
    semilogx(pA, ET, styles{i}, 'LineWidth', 1.8); hold on;
end
ETinf = 1 ./ (1 - p2);
semilogx(pA, ETinf, 'k--', 'LineWidth', 1.2);
grid on; xlabel('p'); ylabel('E[T]');
legend('i = 1','i = 2','i = 3','i \rightarrow \infty','Location','northwest');
title('Expected number of transmissions');

% Numerical check: P(UE | accepted) = p1/(p0+p1) for every i
fprintf('\n=== PART 6: conditional P(UE | accept) ===\n');
for i = [1 2 3 10]
    Pcor = p0 .* geomsum(p2, i);
    Pue  = p1 .* geomsum(p2, i);
    cond = Pue ./ (Pcor + Pue);
    fprintf('i=%2d: max |cond - p1/(p0+p1)| = %.2e\n', i, max(abs(cond - p1./(p0+p1))));
end

%% ========================================================================
%  LOCAL FUNCTIONS
%  ========================================================================
function [C, w, A, dmin, Amin] = codebook(G)
    % All codewords, Hamming weights, spectrum A_0..A_n, dmin and Amin
    k = size(G,1);  n = size(G,2);
    V = dec2bin(0:2^k-1) - '0';          % 2^k x k information vectors
    C = mod(V*G, 2);                     % codebook
    w = sum(C, 2);                       % Hamming weights
    A = accumarray(w+1, 1, [n+1 1]).';   % A(i+1) = A_i
    dmin = find(A(2:end) > 0, 1);        % smallest i >= 1 with A_i > 0
    Amin = A(dmin+1);
end

function pue = PUE(A, p)
    % P_UE(p) = sum_{i=1}^{n} A_i p^i (1-p)^(n-i)
    n = numel(A) - 1;
    pue = zeros(size(p));
    for i = 1:n
        pue = pue + A(i+1) .* p.^i .* (1-p).^(n-i);
    end
end

function ok = verifyH(HT, C, n)
    % Zero syndrome <=> codeword, checked exhaustively over all 2^n vectors
    Y  = dec2bin(0:2^n-1) - '0';
    S  = mod(Y*HT, 2);
    z  = all(S == 0, 2);
    ok = isequal(sortrows(Y(z,:)), sortrows(C));
end

function [nUE, N] = mc_until(G, HT, p, target, batch)
    % Simulate BSC transmissions until `target` undetected errors occur.
    % Returns the number of undetected errors and the exact number N of
    % transmitted codewords at the moment the target is reached.
    k = size(G,1);  n = size(G,2);
    nUE = 0;  N = 0;
    while nUE < target
        V = randi([0 1], batch, k);                  % random info vectors
        C = mod(V*G, 2);                             % encoding
        E = rand(batch, n) < p;                      % BSC error patterns
        Y = xor(C, E);                               % received words
        S = mod(double(Y)*HT, 2);                    % syndromes
        und = all(S == 0, 2) & any(E, 2);            % s = 0 and y ~= c
        cs  = cumsum(und);
        need = target - nUE;
        if cs(end) >= need
            N   = N + find(cs >= need, 1);           % stop exactly at target
            nUE = target;
        else
            N   = N + batch;
            nUE = nUE + cs(end);
        end
    end
end

function S = geomsum(q, i)
    % sum_{j=0}^{i-1} q^j  (computed term by term, safe for q -> 1)
    S = zeros(size(q));
    for j = 0:i-1
        S = S + q.^j;
    end
end
