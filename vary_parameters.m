% parameters to vary
% number of nodes (spatial)
% number of moments
% integration tolerance

% qdot_lw
% coalescence
% nucleation rate
% nucleation distribution


% first, number of nodes:

N_nodes = [10 20 50 100 200 500];
N_mom = [4 6 8 10];
rel_tol = 10.^[-1 -2 -3 -4 -5 -6];

C_qdot_lw = [0.5 1 2 4 8]*1e-4;
C_coalescence = logspace(-5,0,6);
C_coalescence = [C_coalescence(:) ones(length(C_coalescence),2)];
C_nuc_rate = [ 1.5 3 6 12 24]*1e4;

inputs.rel_tol = 1e-4;
inputs.N_mom = 6;
inputs.N_nodes = 50;

inputs.C_qdot_lw = 2e-4;
inputs.C_coalescence = [0 5e-6 5e-2]; % collision efficiency, laminar shear, turbulence
inputs.C_nuc_rate = 6e4;

for i = 1:length(C_nuc_rate);
    inputs.C_nuc_rate = C_nuc_rate(i);
    try
        bubble_growth_1D(inputs)
    end
end

% for i = 1:length(C_coalescence(:,1));
%     inputs.C_coalescence = C_coalescence(i,:);
%     try
%         bubble_growth_1D(inputs)
%     end
% end

% for i = 1:length(C_qdot_lw);
%     inputs.C_qdot_lw = C_qdot_lw(i);
%     try
%         bubble_growth_1D(inputs)
%     end
% end


% bubble_growth_1D(inputs);

% 
% for i = 1:length(N_nodes)
%     inputs.N_nodes = N_nodes(i);
%     try
%     bubble_growth_1D(inputs);
%     end
% end
% 
% inputs.rel_tol = 1e-4;
% inputs.N_nodes = 50;
% 
% for i = 1:length(N_mom)
%     inputs.N_mom = N_mom(i);
%     try
%     bubble_growth_1D(inputs);
%     end
% end
% 
% inputs.N_nodes = 50;
% inputs.N_mom = 6;
% 
% for i = 1:length(rel_tol)
%     inputs.rel_tol = rel_tol(i);
%     try
%     bubble_growth_1D(inputs);
%     end
% end