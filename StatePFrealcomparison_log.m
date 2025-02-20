clc; clear; close all;

% Define n values to test
n_values = 1:15;

% The Toeplitz matrix parameters
r = 2; % Adjust r as needed
T = -r:r;
t_values = ones(2*r+1,1); % Assign all t_k = 1

% Compute T' = T ∩ [-n+1, n-1]
T_prime = T; % Since all valid t_k are 1, we can keep the full range

% Compute new bandwidth r' from T'
r_prime = (length(T_prime) - 1) / 2;

% Compute #states
num_states = nchoosek(2 * r_prime, r_prime); % \binom{2r}{r}

% Construct the state space V
V = cell(num_states, 1);
state_index = 1;

% Generate all states in the correct order, starting from 1100 when r'=2
all_states = nchoosek(1:(2*r_prime), r_prime);

for i = 1:size(all_states,1)
    state = zeros(1, 2*r_prime);
    state(all_states(i, :)) = 1;
    V{state_index} = state;
    state_index = state_index + 1;
end

V = sortrows(cell2mat(V), 'descend');
V = mat2cell(V, ones(size(V,1),1), size(V,2));

% Construct the adjacency matrix A(G_T)
A = zeros(num_states, num_states);

for i = 1:num_states
    original_state = V{i};
    extended_state = [original_state, 0]; % b0
    
    for pos = 1:(2*r_prime+1)
        new_state = extended_state;
        
        if new_state(pos) == 0
            new_state(pos) = 1; % b0+e_k, pos is every try of k
            shifted_state = new_state(2:end); % L(bo+e_k)
            
            % Check if shifted_state matches any valid state in V
            for j = 1:num_states
                if isequal(shifted_state, V{j})
                    k = pos - (r_prime + 1);
                    
                    if ismember(k, T_prime)
                        A(i, j) = 1; % Assign weight 1
                    end
                end
            end
        end
    end
end

% Compute spectral radius (largest eigenvalue)
lambda_max = max(abs(eig(A)));
lambda_values = lambda_max .^ n_values;

perm_approx_values = zeros(size(n_values));
perm_exact_values = zeros(size(n_values));

for idx = 1:length(n_values)
    n = n_values(idx);
    
    % Construct the Toeplitz matrix for each n
    ToeplitzA = zeros(n, n);
    for i = 1:n
        for j = 1:n
            if (j - i) >= -r_prime && (j - i) <= r_prime
                ToeplitzA(i, j) = 1; % All t_k are 1
            end
        end
    end
    
    % Compute matrix power A^n and get A^n(1,1)
    An = A^n;
    perm_A_approx = An(1, 1);
    perm_approx_values(idx) = perm_A_approx;
    
    fprintf('n=%d, perm_approx_values done; ', n);
    drawnow;
    
    % Compute real permanent
    perm_A_exact = perm(ToeplitzA, n);
    perm_exact_values(idx) = perm_A_exact;
    
    fprintf('perm_exact_values done; lambda_values done\n');
    drawnow;
end

% Plot the results
figure;
hold on;
plot(n_values, perm_exact_values, 'rx-', 'LineWidth', 2, 'MarkerSize', 8); % Exact Permanent
plot(n_values, perm_approx_values, 'b^--', 'LineWidth', 2, 'MarkerSize', 8); % Approximate Permanent
plot(n_values, lambda_values, 'gs-', 'LineWidth', 2, 'MarkerSize', 8); % Spectral Radius \lambda^n

% Set log scale for Y-axis
set(gca, 'YScale', 'log');

xlabel('Matrix Size n', 'FontSize', 14);
ylabel('Permanent Values (log scale)', 'FontSize', 14);
legend({'Exact Permanent', 'Approximate Permanent (A^n(1,1))', 'Spectral Radius \lambda^n'}, ...
    'Location', 'NorthWest', 'FontSize', 12);
title(sprintf('Comparison of Permanent Computation Methods (r=%d)', r), 'FontSize', 14);
grid on;

% Define the results folder path
results_folder = '/Users/binghong/Documents/MATLAB/Toeplitz-Matrix-Adjacency-Eigen/results';

% Ensure the folder exists
if ~exist(results_folder, 'dir')
    mkdir(results_folder);
end

% Define the filename and save the figure
filename = fullfile(results_folder, sprintf('perm_comparison_r%d_n%d-%d_log.fig', r, min(n_values), max(n_values)));
savefig(filename);

disp(['Plot saved as: ', filename]);
