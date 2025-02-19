clc; clear; close all;

% The Toeplitz matrix parameters
n = 6;
r = 2;
T = -r:r;
t_values = T;
% t_values = ones(2*r+1,1); % Assign all t_k = 1

% Compute T'=T \cap [-n+1, n-1])
T_prime = T(T >= -n+1 & T <= n-1);

% New bandwidth r' from T'
r_prime = (length(T_prime) - 1) / 2;

% Compute #states
num_states = nchoosek(2 * r_prime, r_prime); % \binom{2r}{r}

% Construct the Toeplitz matrix
ToeplitzA = zeros(n, n);

for i = 1:n
    for j = 1:n
        if (j - i) >= -r_prime && (j - i) <= r_prime
            ToeplitzA(i, j) = t_values(j - i + r + 1);
        end
    end
end

disp('Toeplitz matrix A:');
disp(ToeplitzA);

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
                        t_k_index = find(T == k); % Find the index in T
                        t_k = t_values(t_k_index); % Retrieve correct t_k
                        A(i, j) = t_k;
                    end
                end
            end
        end
    end
end

% Compute matrix power A(G_T)^n and get A(G_T)^n(1,1)
An = A^n;
perm_A_approx = An(1, 1);

% Compute spectral radius (largest eigenvalue)
lambda_max = max(abs(eig(A)));

% Compute real permanent using self-built standard perm(A,n)
perm_A_exact = perm(ToeplitzA, n);

% Display results
disp('Toeplitz matrix:');
disp(ToeplitzA);
disp('Adjacency matrix A(G_T):');
disp(A);
disp('A(G_T)^n:');
disp(An);
disp(['Per(A) (A^n(1,1)): ', num2str(perm_A_approx)]);
disp(['Real Per(A): ', num2str(perm_A_exact)]);
disp(['Spectral radius (growth rate): ', num2str(lambda_max)]);
