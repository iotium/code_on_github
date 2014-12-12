% evaluate effects of C_rdot and C_nuc_rate

A_inj = 7e-7;
specified_case = 6;
C_rdot = 1;
C_nuc_rate = 1;
n_nuc_freq = 3;

C_death_rate = logspace(11,20,9);
n_death_rate = linspace(2,5,3);
alpha_lim = linspace(0.05,1,5);

m = 0;

for i = 1:length(C_death_rate)
    for j = 1:length(n_death_rate)
        for k = 1:length(alpha_lim)
            m = m + 1;
        
%     A_inj = varargin{1};
%     specified_case = varargin{2};
%     E = 1;
%     constants.C_rdot = varargin{3};
%     constants.C_nuc_rate = varargin{4};
%     constants.n_nuc_freq = varargin{5};
%     constants.checking_gauss_error = 0;
%     constants.C_death_rate = varargin{6};
%     constants.n_death_rate = varargin{7};
%     constants.alpha_lim = varargin{8};

    [t{m}, P{m}] = bubble_growth(A_inj, specified_case, C_rdot, C_nuc_rate, n_nuc_freq, C_death_rate(i), n_death_rate(j), alpha_lim(k));
    save eval_death_rate_data.mat
    end
    end
end

