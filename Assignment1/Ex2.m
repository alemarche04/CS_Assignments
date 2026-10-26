clear; clc; close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1. Choose a primitive binary polynomial p(D) of degree r
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% p(D) = D^3 + D + 1
% r = 3


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2. Construct the generator polynomial g(D) = (D+1)p(D)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% g(D) = (D + 1)(D^3 + D + 1) = D^4 + D^3 + D^2 + 1
r = 3;
g = 'z^4 + z^3 + z^2 + 1';
crc_encoder  = comm.CRCGenerator('Polynomial', g);
crc_detector = comm.CRCDetector('Polynomial', g);

k = 2; n = k + r + 1; m = 2^k;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3. Generate all 2^k information vectors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dec2bin - Convert decimal integer to its binary representation
V = dec2bin(0:2^k-1) - '0';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 4. Encode all sequences systematically using g(D)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% obtaining codewords of length n = k + r + 1
Codebook = zeros(m, n);
for i = 1:m
    reset(crc_encoder);
    c = crc_encoder(V(i, :).');
    Codebook(i, :) = c.';
end

% Verify the Codebook with the CRC detector
for i = 1:m
    reset(crc_detector);
    [~, e] = crc_detector(Codebook(i, :).');
    assert(e == 0, 'Codeword %d rejected', i);
end
fprintf('All codewords passed the CRC control\n');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 5. Construct a table 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% containing the four information vectors v, the corresponding information 
% polynomials v(D), the codeword polynomials c(D), the codeword vectors c, 
% and their Hamming weights

w_H = sum(Codebook, 2); % Hamming weights

vD = strings(m, 1);  
cD = strings(m, 1);

for i = 1:m
    vD(i) = poly_str(V(i, :));
    cD(i) = poly_str(Codebook(i, :));
end

T = table(V, vD, Codebook, cD, w_H, ...
    'VariableNames', {'v', 'v(D)','c', 'c(D)', 'Hamming Weights'});
disp(T);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 6. Compute the weight distribution coefficients Ai
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A = zeros(1, n + 1);
for i = 0:n
    A(i + 1) = sum(w_H == i);
end
fprintf('Weight distribution:\n');
for i = 0:n
    fprintf('  A_%d = %d\n', i, A(i + 1));
end
assert(sum(A) == m, 'Sum of A_i must be 2^k');

% plot weight distribution
figure;
stem(0:n, A, 'filled', 'LineWidth', 1.5);
xlabel('Hamming weight i', Interpreter='latex'); 
ylabel('$A_i$', Interpreter='latex');
set(gca, Fontsize = 20);
grid on;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 7. Comment on the results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% particularly on the minimum Hamming weight and the possible codeword 
% weights

dmin = min(w_H(w_H > 0));
fprintf('\nd_min = %d\n', dmin);
fprintf('A_i = %s\n', mat2str(find(A > 0) - 1));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 8. Increase k by one and repeat the procedure 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% until codewords of Hamming weight 2 appear

k = 2;
Atab = {};
while true
    n = k + r + 1;  m = 2^k;
    V = dec2bin(0:m-1) - '0';

    crc_encoder  = comm.CRCGenerator('Polynomial', g);
    crc_detector = comm.CRCDetector('Polynomial', g);

    Codebook = zeros(m, n);
    for i = 1:m
        reset(crc_encoder);
        c = crc_encoder(V(i, :).');
        Codebook(i, :) = c.';
    end

    w_H = sum(Codebook, 2);
    
    A = zeros(1, n + 1);
    for i = 0:n
        A(i + 1) = sum(w_H == i);
    end
    assert(sum(A) == m);

    Atab{end+1} = struct('k', k, 'n', n, 'A', A);   %#ok<SAGROW>
    fprintf('\nk = %2d, n = %2d : d_min = %d', k, n, min(w_H(w_H>0)));

    if A(3) > 0
        fprintf('\n\nCodewords with Hamming weight 2 for k = %d\n\n', k);
        break;
    end
    k = k + 1;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 9. Present a table with weight distribution coefficients for each k
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nmax = Atab{end}.n;
Amat = zeros(numel(Atab), nmax + 1);
for t = 1:numel(Atab)
    a = Atab{t}.A;
    Amat(t, 1:numel(a)) = a;
end
names = compose('A_%d', 0:nmax);
T = array2table(Amat, 'VariableNames', names);
T = addvars(T, cellfun(@(s) s.k, Atab).', cellfun(@(s) s.n, Atab).', ...
    'Before', 1, 'NewVariableNames', {'k', 'n'});
disp(T);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 16. Verify results obtained for points 14 and 15
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k = 4; % Information length

% Generator polynomial g(D) = D^4 + D^3 + D^2 + 1
g_coeffs    = [1 1 1 0 1]; 
r           = length(g_coeffs) - 1;     % deg(g(D)) = 4
n           = k + r;                    % Codeword length = 8

v = [1; 1; 0; 1]; % information vector

% CRC encoding
crc_encoder     = comm.CRCGenerator('Polynomial', 'z^4 + z^3 + z^2 + 1');
crc_detector    = comm.CRCDetector('Polynomial', 'z^4 + z^3 + z^2 + 1');

codeword        = crc_encoder(v);
fprintf('Transmitted Codeword c = [ %s ]\n', num2str(codeword'));

% apply error patterns
% pattern e1: burst of length 3 (detected)
e1 = [0; 0; 1; 0; 1; 0; 0; 0]; 

% pattern e2: burst equal to g(D) (undetected)
e2 = [0; 0; 0; 1; 1; 1; 0; 1]; 

y1 = mod(codeword + e1, 2); % received vector 1
y2 = mod(codeword + e2, 2); % received vector 2

% error detection verification
[~, err1] = crc_detector(y1);
[~, err2] = crc_detector(y2);

fprintf('\n--- ERROR DETECTION RESULTS ---\n');
fprintf('Error Pattern e1 (Burst L = 3):\n');
fprintf('Received vector y1 = [ %s ]\n', num2str(y1'));
fprintf('Detector flag err1 = %d => ', err1);
if err1 == 1
    fprintf('Error detected\n');
else
    fprintf('Undetected error\n');
end

fprintf('\nError Pattern e2 (Burst = g(D)):\n');
fprintf('Received vector y2 = [ %s ]\n', num2str(y2'));
fprintf('Detector flag err2 = %d => ', err2);
if err2 == 1
    fprintf('Error detected\n');
else
    fprintf('Undetected error\n');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% auxiliary function
function s = poly_str(c)
    L  = numel(c); % number of elements in c
    pw = (L-1:-1:0); % decreasing order
    pw = pw(c(:).' == 1);
    if isempty(pw)
        s = "0"; 
        return; 
    end

    t = strings(1, numel(pw));

    for j = 1:numel(pw)

        if pw(j) == 0
            t(j) = "1";
        elseif pw(j) == 1
            t(j) = "D";
        else
            t(j) = "D^" + pw(j);
        end

    end

    s = strjoin(t, " + ");
end