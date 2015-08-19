% parameters to vary
% number of nodes (spatial)
% number of moments
% integration tolerance

% first, number of nodes:

N_nodes = [10 20 50 100 200 500];
N_mom = [4 6 8 10];
rel_tol = 10.^[-1 -2 -3 -4 -5 -6];


inputs.rel_tol = 1e-4;
inputs.N_mom = 6;
inputs.N_nodes = 50;

bubble_growth_1D(inputs);

for i = 1:length(N_nodes)
    inputs.N_nodes = N_nodes(1);
    try
    bubble_growth_1D(inputs);
    end
end

inputs.rel_tol = 1e-4;
inputs.N_nodes = 50;

for i = 1:length(N_mom)
    inputs.N_mom = N_mom(i);
    try
    bubble_growth_1D(inputs);
    end
end

inputs.N_nodes = 50;
inputs.N_mom = 6;

for i = 1:length(rel_tol)
    inputs.rel_tol = rel_tol(i);
    try
    bubble_growth_1D(inputs);
    end
end