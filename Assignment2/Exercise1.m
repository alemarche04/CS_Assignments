clc; clear; close all
rng(22);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

m = 7; % number of cells
L = 2^m - 1;

% primitive polynomials
p1 = [7 4 0];
p2 = [7 1 0];
initCond = [0 0 0 0 0 0 1];


% generate gold sequence with index i >= 0
goldseq = comm.GoldSequence('FirstPolynomial',p1, ...
                            'SecondPolynomial', p2, ...
                            'FirstInitialConditions', initCond, ...
                            'SecondInitialConditions', initCond, ...
                            'Index', 1, ... % index i = 1
                            'SamplesPerFrame', L);

x1  = goldseq()';
x1b = 1 - 2*x1; % bipolar sequence

% cyclic autocorrelation
R1 = zeros(1, L);

for tau = 0:L-1
    R1(tau+1) = sum(x1b .* circshift(x1b, -tau)); 
end

% plot
figure
plot(0:L-1, R1, '.-', MarkerSize=15, LineWidth=1.5)
xlabel('$\tau$', 'Interpreter','latex')
ylabel('Autocorrelation', 'Interpreter','latex')
title('\textbf{Autocorrelation of Gold sequence (index = 1)}', 'Interpreter','latex')
xlim([-1 L]); grid on
set(gca, fontsize = 18)


% generate gold sequence with index i < 0
goldseq = comm.GoldSequence('FirstPolynomial',p1, ...
    'SecondPolynomial', p2, ...
    'FirstInitialConditions', initCond, ...
    'SecondInitialConditions', initCond, ...
    'Index', -1, ... % index i = -1 => c2
    'SamplesPerFrame', L);

x1_neg  = goldseq()';
x1b_neg = 1 - 2*x1_neg; % bipolar sequence

% cyclic autocorrelation
R1_neg = zeros(1, L);

for tau = 0:L-1
    R1_neg(tau+1) = sum(x1b_neg .* circshift(x1b_neg, -tau)); 
end

% plot
figure;
plot(0:L-1, R1_neg, '.-', MarkerSize=15, LineWidth=1.5)
xlabel('$\tau$', 'Interpreter','latex')
ylabel('Autocorrelation', 'Interpreter','latex')
title('\textbf{Autocorrelation of Gold sequence (index = -1)}', 'Interpreter','latex')
xlim([-1 L]); grid on
set(gca, fontsize = 18)


goldseq = comm.GoldSequence('FirstPolynomial',p1, ...
    'SecondPolynomial', p2, ...
    'FirstInitialConditions', initCond, ...
    'SecondInitialConditions', initCond, ...
    'Index', -2, ... % index i = -2 => c1
    'SamplesPerFrame', L);

x2_neg  = goldseq()';
x2b_neg = 1 - 2*x2_neg; % bipolar sequence

% cyclic autocorrelation
R2_neg = zeros(1, L);

for tau = 0:L-1
    R2_neg(tau+1) = sum(x2b_neg .* circshift(x2b_neg, -tau)); 
end

% plot
figure;
plot(0:L-1, R2_neg, '.-', MarkerSize=15, LineWidth=1.5)
xlabel('$\tau$', 'Interpreter','latex')
ylabel('Autocorrelation', 'Interpreter','latex')
title('\textbf{Autocorrelation of Gold sequence (index = -2)}', 'Interpreter','latex')
xlim([-1 L]); grid on
set(gca, fontsize = 18)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

indexes = randperm(L, 2) - 1;

idx1 = indexes(1);
idx2 = indexes(2);

fprintf('Gold sequence indexes: %d and %d\n', idx1, idx2)

% generate the first gold sequence
gold1 = comm.GoldSequence('FirstPolynomial',p1, ...
    'SecondPolynomial', p2, ...
    'FirstInitialConditions', initCond, ...
    'SecondInitialConditions', initCond, ...
    'Index', idx1, ...
    'SamplesPerFrame', L);

x1 = gold1()';
x1b = 1 - 2*x1;


% generate the second gold sequence
gold2 = comm.GoldSequence('FirstPolynomial',p1, ...
    'SecondPolynomial', p2, ...
    'FirstInitialConditions', initCond, ...
    'SecondInitialConditions', initCond, ...
    'Index', idx2, ...
    'SamplesPerFrame', L);

x2 = gold2()';
x2b = 1 - 2*x2;


% cyclic cross-correlation
CC = zeros(1, L);

for tau = 0:L-1
    CC(tau+1) = sum(x1b .* circshift(x2b, -tau)); 
end

% plot
figure;
plot(0:L-1, CC, '.-', MarkerSize=15, LineWidth=1.5)
xlabel('$\tau$', 'Interpreter','latex')
ylabel('Cross-correlation', 'Interpreter','latex')
title('\textbf{Cross-correlation of two Gold sequences}', 'Interpreter','latex')
xlim([-1 L]); ylim([-20 20]); grid on
set(gca, fontsize = 18)


% primitive polynomials that do not form a preferred pair
p1 = [7 4 0];
p2 = [7 3 0];

% first m-sequence
pn1 = comm.PNSequence( ...
    'Polynomial', p1, ...
    'InitialConditions', initCond, ...
    'SamplesPerFrame', L);

u1  = pn1()';
u1b = 1 - 2*u1;


% second m-sequence
pn2 = comm.PNSequence( ...
    'Polynomial', p2, ...
    'InitialConditions', initCond, ...
    'SamplesPerFrame', L);

u2  = pn2()';
u2b = 1 - 2*u2;


% cyclic cross-correlation
CC = zeros(1, L);

for tau = 0:L-1
    CC(tau+1) = sum(u1b .* circshift(u2b, -tau)); 
end

% plot
figure;
plot(0:L-1, CC, '.-', MarkerSize=15, LineWidth=1.5)
xlabel('$\tau$', 'Interpreter','latex')
ylabel('Cross-correlation', 'Interpreter','latex')
title('\textbf{Cross-correlation of two non-preferred m-sequences}', 'Interpreter','latex')
xlim([-1 L]); grid on
set(gca, fontsize = 18)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part 3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% calculate the Welch bound
Kvalues = 1:200;

Welch_general = zeros(size(Kvalues));
Welch_simple = zeros(size(Kvalues));
Welch_approx = zeros(size(Kvalues));

for k = 1:length(Kvalues)

    Kv = Kvalues(k);

    % general Welch bound for s = 1,...,10
    values = zeros(1,10);

    for s = 1:10

        argument = (1/(Kv*L-1)) * ...
            (Kv*L/nchoosek(L+s-1,s) - 1);

        argument = max(0,argument);

        values(s) = L*argument^(1/(2*s));

    end

    Welch_general(k) = ceil(max(values));

    % simplified Welch bound (s = 1)
    Welch_simple(k) = ceil( ...
        L*sqrt((Kv-1)/(Kv*L-1)));

    % large K approximation
    Welch_approx(k) = ceil(sqrt(L));

end

% calculate the Sidelnikov bound

Sidelnikov = zeros(size(Kvalues));

for k = 1:length(Kvalues)

    Kv = Kvalues(k);

    values = zeros(1,11);

    for s = 0:10

        argument = (2*s+1)*(L-s) + ...
            s*(s+1)/2 - ...
            (2^s * L^(2*s+1)) / ...
            (Kv*factorial(2*s)*nchoosek(L,s));

        argument = max(0,argument);

        values(s+1) = sqrt(argument);

    end

    Sidelnikov(k) = ceil(max(values));

end

figure;
hold on;
plot(Kvalues,Sidelnikov,'r','LineWidth',1.5);
plot(Kvalues,Welch_general,'b-o','LineWidth',1.5);
plot(Kvalues,Welch_simple,'c','LineWidth',1.5);
plot(Kvalues,Welch_approx,'k','LineWidth',1.5);
hold off; 
legend('Sidelnikov binary', ...
    'Welch original', ...
    'Welch simplified formula', ...
    'Welch approximation', ...
    'Interpreter','latex', ...
    'Location', 'southeast');
grid on; ylim([0 20]);
title('\textbf{Welch and Sidelnikov bounds (N = 127)}', 'Interpreter','latex')
set(gca, fontsize = 18)



% compute rM considering the entire set of K = 129 Gold sequences
K = 129;

Cb = zeros(K,L);

% Gold sequence indexes: -2, -1, 0, ..., 126
indexes = [-2 -1 0:L-1];

p1 = [7 4 0];
p2 = [7 1 0];
initCond = [0 0 0 0 0 0 1];

for k = 1:K

    goldseq = comm.GoldSequence( ...
        'FirstPolynomial', p1, ...
        'SecondPolynomial', p2, ...
        'FirstInitialConditions', initCond, ...
        'SecondInitialConditions', initCond, ...
        'Index', indexes(k), ...
        'SamplesPerFrame', L);

    x = goldseq()';

    % convert binary to bipolar
    Cb(k,:) = 1 - 2*x;

end

% compute rA (maximal out-of-phase cyclic auto-correlation magnitude) 
% and rC (maximal cyclic cross-correlation magnitude)

rA = 0;
rC = 0;

for tau = 0:L-1

    % cyclically shift all sequences
    Cshift = circshift(Cb,[0 -tau]);

    % cross-correlation matrix
    R = Cb * Cshift';

    % cross-correlation: different sequences
    Rcross = abs(R);
    Rcross(1:K+1:end) = 0;

    rC = max(rC,max(Rcross(:)));

    % autocorrelation: same sequence, nonzero shifts
    if tau > 0
        rA = max(rA,max(abs(diag(R))));
    end

end

rM = max(rA,rC); % maximal cyclic correlation magnitude

fprintf('\n_____ Gold family correlation statistics _____\n');
fprintf('rA = %d\n',rA);
fprintf('rC = %d\n',rC);
fprintf('rM = %d\n',rM);


% Add actual Gold family result
figure;
hold on;
plot(Kvalues,Sidelnikov,'r','LineWidth',1.5);
plot(Kvalues,Welch_general,'b-o','LineWidth',1.5);
plot(Kvalues,Welch_simple,'c','LineWidth',1.5);
plot(Kvalues,Welch_approx,'k','LineWidth',1.5);
plot(K,rM,'ko','MarkerFaceColor','k', 'MarkerSize',8);
hold off; 
legend('Sidelnikov binary', ...
    'Welch original', ...
    'Welch simplified formula', ...
    'Welch approximation', ...
    'Gold family $r_M$', ...
    'Interpreter','latex', ...
    'Location', 'southeast');
grid on; ylim([0 20]);
title('\textbf{Welch and Sidelnikov bounds (N = 127)}', 'Interpreter','latex')
set(gca, fontsize = 18)


% Display results at K = 129

fprintf('\n_____ Bounds at K = 129 _____\n');

fprintf('Welch general = %d\n',Welch_general(129));
fprintf('Welch simplified = %d\n',Welch_simple(129));
fprintf('Welch approximation = %d\n',Welch_approx(129));
fprintf('Sidelnikov = %d\n',Sidelnikov(129));
fprintf('Actual Gold rM = %d\n',rM);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part 4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

m = 7;
L = 2^m - 1;
N = 126;
K = 2^m + 1;

% generate the complete Gold family
Cb = zeros(K,L);

indexes = [-2 -1 0:L-1];

p1 = [7 4 0];
p2 = [7 1 0];
initCond = [0 0 0 0 0 0 1];

for k = 1:K

    goldseq = comm.GoldSequence( ...
        'FirstPolynomial', p1, ...
        'SecondPolynomial', p2, ...
        'FirstInitialConditions', initCond, ...
        'SecondInitialConditions', initCond, ...
        'Index', indexes(k), ...
        'SamplesPerFrame', L);

    x = goldseq()';

    % convert binary to bipolar
    Cb(k,:) = 1 - 2*x;

end

% truncate sequences to length 126
C126 = Cb(:, 1:N);

% compute rA (maximal out-of-phase cyclic auto-correlation magnitude) 
% and rC (maximal cyclic cross-correlation magnitude)

rA = 0;
rC = 0;

for tau = 0:N-1

    % cyclically shift all sequences
    Cshift = circshift(C126,[0 -tau]);

    % cross-correlation matrix
    R = C126 * Cshift';

    % cross-correlation: different sequences
    Rcross = abs(R);
    Rcross(1:K+1:end) = 0;

    rC = max(rC,max(Rcross(:)));

    % autocorrelation excluding zero shift
    if tau > 0
        rA = max(rA,max(abs(diag(R))));
    end

end

rM = max(rA,rC); % maximal cyclic correlation magnitude

fprintf('\n_____ Truncated Gold family correlation statistics _____\n');
fprintf('rA = %d\n',rA);
fprintf('rC = %d\n',rC);
fprintf('rM = %d\n',rM);


% calculate the Welch bound
Kvalues = 1:200;

Welch_general   = zeros(size(Kvalues));
Welch_simple    = zeros(size(Kvalues));
Welch_approx    = zeros(size(Kvalues));

for k = 1:length(Kvalues)

    Kv = Kvalues(k);

    % general Welch bound for s = 1,...,10
    values = zeros(1,10);

    for s = 1:10

        argument = (1/(Kv*N-1)) * ...
            (Kv*N/nchoosek(N+s-1,s) - 1);

        argument = max(0,argument);

        values(s) = N*argument^(1/(2*s));

    end

    Welch_general(k) = ceil(max(values));

    % simplified Welch bound (s = 1)
    Welch_simple(k) = ceil( ...
        N*sqrt((Kv-1)/(Kv*N-1)));

    % large K approximation
    Welch_approx(k) = ceil(sqrt(N));

end

% calculate the Sidelnikov bound

Sidelnikov = zeros(size(Kvalues));

for k = 1:length(Kvalues)

    Kv = Kvalues(k);

    values = zeros(1,11);

    for s = 0:10

        argument = (2*s+1)*(N-s) + ...
            s*(s+1)/2 - ...
            (2^s * N^(2*s+1)) / ...
            (Kv*factorial(2*s)*nchoosek(N,s));

        argument = max(0,argument);

        values(s+1) = sqrt(argument);

    end

    Sidelnikov(k) = ceil(max(values));

end

figure;
hold on;
plot(Kvalues,Sidelnikov,'r','LineWidth',1.5);
plot(Kvalues,Welch_general,'b-o','LineWidth',1.5);
plot(Kvalues,Welch_simple,'c','LineWidth',1.5);
plot(Kvalues,Welch_approx,'k','LineWidth',1.5);
hold off; 
legend('Sidelnikov binary', ...
    'Welch original', ...
    'Welch simplified formula', ...
    'Welch approximation', ...
    'Interpreter','latex', ...
    'Location', 'southeast');
grid on; ylim([0 20]);
title('\textbf{Welch and Sidelnikov bounds (N = 126)}', 'Interpreter','latex')
set(gca, fontsize = 18)

fprintf('\n_____ Bounds at K = 129 _____\n');

fprintf('Welch general = %d\n',Welch_general(129));
fprintf('Welch simplified = %d\n',Welch_simple(129));
fprintf('Welch approximation = %d\n',Welch_approx(129));
fprintf('Sidelnikov = %d\n',Sidelnikov(129));
fprintf('Actual Gold rM = %d\n',rM);


% Add rM to the plot
figure;
hold on;
plot(Kvalues,Sidelnikov,'r','LineWidth',1.5);
plot(Kvalues,Welch_general,'b-o','LineWidth',1.5);
plot(Kvalues,Welch_simple,'c','LineWidth',1.5);
plot(Kvalues,Welch_approx,'k','LineWidth',1.5);
plot(K,rM,'ko','MarkerFaceColor','k', 'MarkerSize',8);
hold off; 
legend('Sidelnikov binary', ...
    'Welch original', ...
    'Welch simplified formula', ...
    'Welch approximation', ...
    'Gold family $r_M$', ...
    'Interpreter','latex', ...
    'Location', 'northwest');
grid on; ylim([0 50]);
title('\textbf{Welch and Sidelnikov bounds (N = 126)}', 'Interpreter','latex')
set(gca, fontsize = 18)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part 5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = 128;
K = 128;

% generate Walsh-Hadamard matrix
H = hadamard(N);

% verify zero-delay orthogonality
R0 = H * H';
expected = N*eye(N); % expected result

orthogonality_error = max(abs(R0(:)-expected(:)));

fprintf('\n_____ Walsh-Hadamard matrix _____\n')
if orthogonality_error == 0
    fprintf('Zero-delay orthogonality satisfied: at zero relative delay, all rows are perfectly orthogonal\n');
else
    fprintf('ERROR: Zero-delay orthogonality not satisfied\n');
end


% compute rA (maximal out-of-phase cyclic auto-correlation magnitude) 
% and rC (maximal cyclic cross-correlation magnitude)
rA = 0;
rC = 0;

for tau = 0:N-1

    % cyclically shift all rows
    Hshift = circshift(H,[0 -tau]);

    % cross-correlation matrix
    R = H * Hshift';

    % cross-correlation between different rows
    Rcross = abs(R);

    % remove diagonal (autocorrelation)
    Rcross(1:K+1:end) = 0;

    rC = max(rC,max(Rcross(:)));

    % autocorrelation excluding zero shift
    if tau > 0
        rA = max(rA,max(abs(diag(R))));
    end

end

rM = max(rA,rC); % maximal cyclic correlation magnitude

fprintf('\n_____ Walsh-Hadamard correlation statistics _____\n');
fprintf('rA = %d\n',rA);
fprintf('rC = %d\n',rC);
fprintf('rM = %d\n',rM);


% calculate the Welch bound
Kvalues = 1:200;

Welch_general   = zeros(size(Kvalues));
Welch_simple    = zeros(size(Kvalues));
Welch_approx    = zeros(size(Kvalues));

for k = 1:length(Kvalues)

    Kv = Kvalues(k);

    % general Welch bound for s = 1,...,10
    values = zeros(1,10);

    for s = 1:10

        argument = (1/(Kv*N-1)) * ...
            (Kv*N/nchoosek(N+s-1,s) - 1);

        argument = max(0,argument);

        values(s) = N*argument^(1/(2*s));

    end

    Welch_general(k) = ceil(max(values));

    % simplified Welch bound (s = 1)
    Welch_simple(k) = ceil( ...
        N*sqrt((Kv-1)/(Kv*N-1)));

    % large K approximation
    Welch_approx(k) = ceil(sqrt(N));

end

% calculate the Sidelnikov bound

Sidelnikov = zeros(size(Kvalues));

for k = 1:length(Kvalues)

    Kv = Kvalues(k);

    values = zeros(1,11);

    for s = 0:10

        argument = (2*s+1)*(N-s) + ...
            s*(s+1)/2 - ...
            (2^s * N^(2*s+1)) / ...
            (Kv*factorial(2*s)*nchoosek(N,s));

        argument = max(0,argument);

        values(s+1) = sqrt(argument);

    end

    Sidelnikov(k) = ceil(max(values));

end


% Add actual Gold family result
figure;
hold on;
plot(Kvalues,Sidelnikov,'r','LineWidth',1.5);
plot(Kvalues,Welch_general,'b-o','LineWidth',1.5);
plot(Kvalues,Welch_simple,'c','LineWidth',1.5);
plot(Kvalues,Welch_approx,'k','LineWidth',1.5);
plot(K,rM,'ko','MarkerFaceColor','k', 'MarkerSize',8);
hold off; 
legend('Sidelnikov binary', ...
    'Welch original', ...
    'Welch simplified formula', ...
    'Welch approximation', ...
    'Gold family $r_M$', ...
    'Interpreter','latex', ...
    'Location', 'northwest');
grid on; ylim([0 130]);
title('\textbf{Welch and Sidelnikov bounds (N = 128)}', 'Interpreter','latex')
set(gca, fontsize = 18)


% Display results at K = 128

fprintf('\n_____ Bounds at K = 129 _____\n');

fprintf('Welch general = %d\n',Welch_general(128));
fprintf('Welch simplified = %d\n',Welch_simple(128));
fprintf('Welch approximation = %d\n',Welch_approx(128));
fprintf('Sidelnikov = %d\n',Sidelnikov(128));
fprintf('Actual Gold rM = %d\n',rM);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Cyclic cross-correlation function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute Welch bounds function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute Sidelnikov bound function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%