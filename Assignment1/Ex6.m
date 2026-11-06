clc
clear
close all

%% Property of m-sequences
% choose a primitive polynomial
primitivePolynomial = [10 9 4 1 0];

% choose a random seed
initialConditions = [1 1 1 0 0 0 1 0 1 0];

L   = 10; % number of cells
Nb  = 2^L-1; % m-sequence length

% sequence of pseudorandom binary numbers using a LFSR
pnSequence = comm.PNSequence('Polynomial', primitivePolynomial, ...
'SamplesPerFrame',Nb,'InitialConditions', initialConditions);

x1 = pnSequence()'; % binary m-sequence

% generate the bipolar sequence of x1
x1_ = 2.*x1 - 1;

%% Sequence
% plot the first 20 bits of the bipolar sequence
figure
% stem(0:19, x1_(1:20), 'filled')
stairs(0:20, x1_(1:21), LineWidth=1.5)
grid on; 
ylim([-2 2]); yticks([-1 0 1]); yticklabels([-1 0 1])
xticks(0:20); xticklabels(0:20)
set(gca, fontsize = 18)

%% Number of 0/1, number of transitions
N1 = sum(x1);
N0 = sum(x1==0);

% consider sequence x1 as periodic
x1_per = [x1 x1(1)];

NT  = sum(diff(x1_per)~=0); % number of transitions
NNT = sum(diff(x1_per) == 0); % number of non-transitions

% write a table with N0, N1, NT and NNT
sequenceStatistics = table(N0, N1, NT, NNT, ...
    'VariableNames', {'N0', 'N1', ...
    'NT', 'NNT'});
disp(sequenceStatistics)

%% Runs
% Count runs of consecutive equal symbols, including the periodic boundary
N       = length(x1);
x1_per  = [x1 x1(1)];

% find transition points between adjacent bits (including periodic wrap-around)
trans = find(diff(x1_per) ~= 0);  % end index of each run

% calculate run lengths across periodic boundaries
runLengths = diff([trans(end) - N, trans]);

% extract the symbol value (0 or 1) corresponding to each run
runValues = x1(trans);

% calculate starting position in x1 for each run (1-based indexing)
trans_prev = [trans(end), trans(1:end-1)];
runStarts  = mod(trans_prev, N) + 1;

% separate Statistics for 1s and 0s
NR              = length(runLengths);
one_runLengths  = runLengths(runValues == 1);
zero_runLengths = runLengths(runValues == 0);

one_runLongest  = max(one_runLengths);
zero_runLongest = max(zero_runLengths);

fprintf('Total Runs (NR)        : %d\n', NR);
fprintf('Total Runs of 1s       : %d\n', sum(numel(one_runLengths)));
fprintf('Total Runs of 0s       : %d\n', sum(numel(zero_runLengths)));
fprintf('Longest Run of 1s      : %d\n', one_runLongest);
fprintf('Longest Run of 0s      : %d\n', zero_runLongest);


% construct Run-Length distribution tables (Lengths 1 to 10)
varNames = {'NR(1)', 'NR(2)', 'NR(3)', 'NR(4)', 'NR(5)', ...
            'NR(6)', 'NR(7)', 'NR(8)', 'NR(9)', 'NR(10)'};

% total runs distribution
NR_counts = histcounts(runLengths, 1:11);
NR_i = array2table(NR_counts, 'VariableNames', varNames);
fprintf('\n\nTotal Runs Distribution N_R(i)\n\n:');
disp(NR_i);

% runs of 1s distribution
varNames1 = {'NR1(1)', 'NR1(2)', 'NR1(3)', 'NR1(4)', 'NR1(5)', ...
             'NR1(6)', 'NR1(7)', 'NR1(8)', 'NR1(9)', 'NR1(10)'};
NR1_counts = histcounts(one_runLengths, 1:11);
NR1_i = array2table(NR1_counts, 'VariableNames', varNames1);
fprintf('\n\nRuns of 1s Distribution N_R1(i):\n\n');
disp(NR1_i);

% runs of 0s distribution
varNames0 = {'NR0(1)', 'NR0(2)', 'NR0(3)', 'NR0(4)', 'NR0(5)', ...
             'NR0(6)', 'NR0(7)', 'NR0(8)', 'NR0(9)', 'NR0(10)'};
NR0_counts = histcounts(zero_runLengths, 1:11);
NR0_i = array2table(NR0_counts, 'VariableNames', varNames0);
fprintf('\n\nRuns of 0s Distribution N_R0(i):\n\n');
disp(NR0_i);

% plot Run Lengths vs Starting Position
figure;
hold on;

is_zero = (runValues == 0);
is_one  = (runValues == 1);

% stem plot for runs of 0s (blue)
stem(runStarts(is_zero), runLengths(is_zero), 'b', 'LineWidth', 1.2, ...
    'Marker', '.', 'MarkerSize', 8, 'DisplayName', 'Runs of 0s');

% stem plot for runs of 1s (magenta)
stem(runStarts(is_one), runLengths(is_one), 'm', 'LineWidth', 1.2, ...
    'Marker', '.', 'MarkerSize', 8, 'DisplayName', 'Runs of 1s');

grid on;
xlabel('Starting position in the sequence');
ylabel('Run length');
legend('Location', 'northeast');
xlim([0, N + 10]);
ylim([0, max(runLengths) + 1]);
set(gca, fontsize=18)
hold off;


%% Autocorrelation
R   = ifft(fft(x1_).*conj(fft(x1_))); % non-normalized periodic autocorr.
R   = fftshift(R);

tau_axis = -(N-1)/2 : (N-1)/2;

figure;
plot(tau_axis, R, 'LineWidth', 1.5);
grid on;
xlabel('\tau');
ylabel('R_{v''}(\tau)');
xlim([-(N-1)/2, (N-1)/2]);
ylim([-50 1100]);
set(gca, fontsize=18)


%% Truncated sequences: Maximum Peak Side Lobe (MPSL)

z_values = 1:200;

MPSL = zeros(size(z_values));

for k = 1:length(z_values)

    z = z_values(k);

    % Remove the last z bits
    xz = x1_(1:N-z);

    % Compute linear (aperiodic) autocorrelation
    R_aperiodic = xcorr(xz);

    % Zero out the main peak at lag tau = 0
    center_idx = length(xz);
    R_aperiodic(center_idx) = 0;

    % Maximum positive side-lobe
    MPSL(k) = max(abs(R_aperiodic));

end

% plot MPSL vs z
figure;
plot(z_values,MPSL,'LineWidth',1.5);
xlabel('z');
ylabel('MPSL(z)');
grid on; axis tight;
set(gca, fontsize=18)


% display some results
selected_z = [1 2 10 23 50 100 150 200];
selected_MPSL = MPSL(selected_z);
MPSLTable = table(selected_z',selected_MPSL', ...
    'VariableNames',{'RemovedBits','MPSL'});
disp(MPSLTable);

%% Truncated sequences: seed
z_fixed = 23;
primitivePolynomial = [10 9 4 1 0];

% seed A: [1 1 1 0 0 0 1 0 1 0];
initialConditions_A = [1 1 1 0 0 0 1 0 1 0];
pnGenA = comm.PNSequence('Polynomial', primitivePolynomial, ...
    'InitialConditions', initialConditions_A, ...
    'SamplesPerFrame', Nb);
vA = pnGenA()';
vA_bipolar = 2*vA - 1;
vA_trunc = vA_bipolar(1:N-z_fixed);

% seed B: [0 0 0 0 0 1 0 1 1 0]
initialConditions_B = [0 0 0 0 0 1 0 1 1 0];
pnGenB = comm.PNSequence('Polynomial', primitivePolynomial, ...
    'InitialConditions', initialConditions_B, ...
    'SamplesPerFrame', Nb);
vB = pnGenB()';
vB_bipolar = 2*vB - 1;
vB_trunc = vB_bipolar(1:N-z_fixed);

% aperiodic autocorrelation
R_trunc_A = xcorr(vA_trunc);
R_trunc_B = xcorr(vB_trunc);
lags = -(length(vA_trunc)-1):(length(vA_trunc)-1);

figure;
subplot(2,1,1);
plot(lags, R_trunc_A, 'b-', 'LineWidth', 1.5);
grid on;
title("$R_{v'_z}(\tau)$ ($z = 23$, Seed A: [1110001010])", 'Interpreter','latex');
xlabel('$\tau$', 'Interpreter','latex'); ylabel("$R_{v'_z}(\tau)$", 'Interpreter','latex');
ylim([-100 1050]);
set(gca, fontsize=22)

subplot(2,1,2);
plot(lags, R_trunc_B, 'm-', 'LineWidth', 1.5);
grid on;
title("$R_{v'_z}(\tau)$ ($z = 23$, Seed B: [0000010110]))", 'Interpreter','latex');
xlabel('$\tau$', 'Interpreter','latex'); ylabel("$R_{v'_z}(\tau)$", 'Interpreter','latex');
ylim([-100 1050]);
set(gca, fontsize=22)
