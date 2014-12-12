% evaluate effects of C_rdot and C_nuc_rate

C_rdot = logspace(-2,2,30);
C_nuc_rate = logspace(-2,2,30);

A_inj = 7e-7;
E = 1;
specified_case = 6;

for i = 1:length(C_rdot)
    for j = 1:length(C_nuc_rate)
        
%         varargout{1} = t(n_min);
% varargout{2} = t(n_peak);
% varargout{3} = P(n_min);
% varargout{4} = P(n_peak);
% 
%     A_inj = varargin{1};
%     specified_case = varargin{2};
%     E = varargin{3};
%     constants.C_rdot = varargin{4};
%     constants.C_nuc_rate = varargin{5};

[t_min(i,j), t_peak(i,j), P_min(i,j), P_peak(i,j)] = ...
    bubble_growth(A_inj, specified_case, E, C_rdot(i), C_nuc_rate(j));

save eval_C_data
    end
end

