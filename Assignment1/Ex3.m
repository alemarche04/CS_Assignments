clc
close all
clear

%% 3b.2
n = 6; k = 3; r = n - k;

% A non-systematic generator matrix is a generator matrix kxn G for a 
% linear block code that does not contain an explicit identity submatrix 
% I_k in a % partitioned form like G = [I_k | P]

% 3x6 generator matrix that is non-systematic and has full rank 3
G = [1, 1, 1, 0, 0, 1;
     0, 1, 1, 1, 0, 0;
     1, 0, 1, 1, 1, 0];

% generate all 2^k = 8 information vectors
% dec2bin - Convert decimal integer to its binary representation
V = dec2bin(0:2^k-1) - '0'; % - '0' => conversion from binary strings to binary values

% encode all information vectors
C = mod(V * G, 2);

% all possible received vectors
Y = dec2bin(0:2^n-1) - '0';

%% 3b.5
% decision regions
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

fprintf('--- Decision Regions D(c) for each codeword ---\n\n');
for i = 1:size(C, 1)
    fprintf('Decision Region D(c_%d) for codeword c_%d = [%s] (%d vectors):\n', ...
        i-1, i-1, num2str(C(i,:)), size(D_c{i}, 1));
    disp(D_c{i});
end


%% 3b.6
% find the all-zero codeword index
zero_idx = find(all(C == 0, 2));
D_zero = D_c{zero_idx};

verified = true; % initialization

fprintf('--- Verification of D(c) = c + D(0) ---\n');
for i = 1:size(C, 1)
    c = C(i, :);

    % compute c + D(0) by adding c (modulo 2) to each vector in D(0)
    c_plus_D0 = mod(D_zero + c, 2);

    % sort both sets to compare them regardless of order
    current_region = sortrows(D_c{i});
    expected_region = sortrows(c_plus_D0);

    is_equal = isequal(current_region, expected_region);
    fprintf('Codeword c_%d: D(c_%d) = c_%d + D(0) -> %s\n', ...
        i-1, i-1, i-1, mat2str(is_equal));

    if ~is_equal
        verified = false;
    end
end

if verified
    fprintf('\nVerification SUCCESSFUL: D(c) = c + D(0) holds for all codewords c in C.\n');
else
    fprintf('\nVerification FAILED for some codewords.\n');
end

%% 3c.1
n = 4; k = 2; r = n - k;

% define nonlinear code C
C = [   0, 0, 0, 0; 
        1, 1, 0, 0;  
        0, 1, 1, 0;  
        1, 1, 1, 1 ]; 

% all possible received vectors
Y = dec2bin(0:2^n-1) - '0';

%% 3c.3
% decision regions
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

fprintf('--- Decision Regions D(c) for each codeword ---\n\n');
for i = 1:size(C, 1)
    fprintf('Decision Region D(c_%d) for codeword c_%d = [%s] (%d vectors):\n', ...
        i-1, i-1, num2str(C(i,:)), size(D_c{i}, 1));
    disp(D_c{i});
end


%% 3c.4
% find the all-zero codeword index
zero_idx = find(all(C == 0, 2));
D_zero = D_c{zero_idx};

verified = true; % initialization

fprintf('--- Verification of D(c) = c + D(0) ---\n');
for i = 1:size(C, 1)
    c = C(i, :);

    % compute c + D(0) by adding c (modulo 2) to each vector in D(0)
    c_plus_D0 = mod(D_zero + c, 2);

    % sort both sets to compare them regardless of order
    current_region = sortrows(D_c{i});
    expected_region = sortrows(c_plus_D0);

    is_equal = isequal(current_region, expected_region);
    fprintf('Codeword c_%d: D(c_%d) = c_%d + D(0) -> %s\n', ...
        i-1, i-1, i-1, mat2str(is_equal));

    if ~is_equal
        verified = false;
    end
end

if verified
    fprintf('\nVerification SUCCESSFUL: D(c) = c + D(0) holds for all codewords c in C.\n');
else
    fprintf('\nVerification FAILED for some codewords.\n');
end