% this is the main file for my tank modeling code
% it can be run by itself (with no arguments)
% or it can be called by another script and arguments are passed to it

% basic outline: 2 functions
% bubble_growth_1D: set up the problem, integrate differential equations,
% plot results
% diffeqns: calculate the d/dt terms, these are passed to bubble_growth_1D
% and integrated

% numerous other functions are called, 
% but only these two are in this m file

function varargout = bubble_growth_1D(varargin)
%--------------------------------------------------------------------------
% set physical parameters of the problem
%--------------------------------------------------------------------------

% feedline parameters
constants.f_feedline = 0.0; % friction factor for the feedline (from the tank to the orifice)
constants.L_feedline = 5 * 0.0254;  % [m] length of the feedline
constants.D_feedline = (.375 - 2*0.049) * 0.0254;   % [m] ID of the feedline


specified_case = 12; % the initial_conditions.m file contains many different 
% sets of initial conditions that correspond to experimental data sets 
% that I compare to. This number specifies which one of these are to be
% used. Open up that m file to see more. If 0 is entered, the data below
% are used.

% initial conditions

% parameters I can play around with
if specified_case == 0
    % N2O test 2 from my data
    Ti = 283.7;           % [K] initial temperature
    fill_level = 0.80;        % [] initial fill_level ratio (by volume)
    %         E = 2.1e4;          % [] heat transfer multiplier
    V_tank = 1.80e-4;   % [m^3] tank volume
    L_tank = 0.356;     % [m] tank length
    Cd = 1;         % [] injector Cd
    %         A_inj = 3.15e-7;       % [m^2] injector area
    Po = 1e5;           % [Pa] external pressure
    T_air = 293;        % [K] air temperature
    rho_w = 1360;       % [kg/m^3] density of wall material (polycarb)
    cv_w = 1250;        % [J/kg.K] specific heat of wall (polycarb)
    t_w = 0.0254*1/4;   % [m] wall thickness
    D_tank = sqrt(4/pi*V_tank/L_tank);
    % [m] tank diameter
    k_w = 0.195;          % [W/m.K] thermal conductivity of wall
else
    % use one of the specified cases in initial_conditions.m
    [Ti, fill_level, V_tank, L_tank, ...
        A_inj, Cd, Po, T_air, rho_w, cv_w, t_w, D_tank, k_w, fluid] = initial_conditions(specified_case);
    
end

%--------------------------------------------------------------------------
% set model parameters
%--------------------------------------------------------------------------

constants.use_numjac = 0;   % whether or not to use the numjac function to 
% numerically calculate the jacobian (0 = use my own code)

% options for fsolve, used to solve a variety of nonlinear equations
constants.fsolve_options = optimset('TolX',1e-12,'Display','off');

% use either REFPROP, or interpolate from a table (PDT)
constants.property_source = 'PDT';
% constants.property_source = 'refprop';

% nucleation model: Shin & Jones, or Alamgir & Lienhard
constants.nuc_model = 'SJ';
% constants.nuc_model = 'AL';

% Jensen & Memmel departure diameter correlation. They came up with one
% that includes the wall superheat, and one that does not
constants.r_dep_expression = 'with superheat';
% constants.r_dep_expression = 'without superheat';

% nucleation density expressions. Hibiki and Ishii is more modern, Shin and
% Jones may be more consistent with other aspects of the nucleation model
% I'm using
constants.nuc_density_expression = 'hibiki and ishii';
% constants.nuc_density_expression = 'shin and jones';

% nucleation frequency expressions. Shin and Jones, or a very loose
% correlation based on the experiments of Saddy and Jameson
% constants.nuc_frequency_expression = 'saddy and jameson';
constants.nuc_frequency_expression = 'shin and jones';

% which norm to use when controlling error.
constants.error_norm = 'L-infinity';
% constants.error_norm = 'L-2';
% constants.error_norm = 'L-1';

constants.include_rdot_rhodot = 0;  % include bubble size changes that occur due to gas density changing
constants.include_hysteresis = 1;   % include hysteresis in the bubble nucleation process
constants.include_u_bulk = 0;       % include a bulk velocity based on the liquid flowing out the orifice
constants.ADQMOM_p = 1;             % set the scale factor for ADQMOM
constants.adaptive_mesh_refinement = 0; % include adaptive mesh refinement

% parameters needed for adaptive mesh refinement
max_spatial_rel_delta_tol = 2;      % max difference between two nodes before refinement occurs
min_spatial_rel_delta_tol = 0.025;  % min difference between two nodes before coarsening occurs
max_L = 0.1;                        % max distance between nodes (relative to tank length)

ADQMOM = 'off'; % whether or not to use the ADQMOM technique 
% (not fully supported anymore, so leave 'off')

save_filename = 'bubble_sim_data';  % file name to save data to

% ODE integration scheme
ode_solver = 'TRBDF2'; % [] options:
% adaptive:
% 'RKF' for runge-kutta-fehlberg
% 'CK' for cash-karp
% 'DP' for dormand-prince
% 'BS' for 2/3rd order bogacki-shampine
% 'ROCK2' for the rock2 orthogonal chebyshev runge kutta
% 'ROS3P' for ROS3P 3rd order rosenbrock
% 'TRBDF2' for trapezoid/backwards difference formula 2

% not adaptive:
% 'euler' for 1st order euler
% 'RK4' for 4th order runge-kutta

ROCK2_stages = 20;  % number of stages in the rock2 orthogonal chebyshev 
% runge kutta scheme (if that was selected above)

constants.sh_max = 1.25; % max relative increase in step size
sh_min = 0.1; % min relative decrease in step size

N_nodes = 32;   % number of nodes in the vertical direction
N_mom = 4;  % number of moments tracked in DQMOM
rel_tol = 1e-3;     % [] max relative error allowed in adaptive scheme
constants.C_qdot_lw = 2e-5; % multiplies the heat transfer from wall to liquid
constants.C_coalescence = [0.1 1 1]; % collision efficiency,
% laminar shear, turbulence (buoyancy is set at 1)
constants.C_nuc_rate = 1e-8; % multiplies the nucleation rate

newton_tol = 0.1; % (tolerance for quasi newton iteration, relative to rel_tol)
max_iter = 5; % max number of quasi newton iterations

constants.phi = 0.1;    % parameter of the alamgir & lienhard nucleation rate expression
constants.CD_churn_turb_boundary = 0.3; % void fraction at which the drag correlation switches to churn-turbulent

N_rw = 25;  % number of nodes in the radial discretization of the wall for conduction

force_constant_h = 0;   % use a constant time step rather than adaptive

h = 1e-5;           % [s] initial time step
abs_tol = 1e9;     % [] max absolute error allowed in adaptive scheme
h_min = 1e-14;      % [s] min allowable time step
t_end = 1e3;         % [s] end time (if LRO doesn't happen first)
LRO_tol = 2e-2;     % [s] tolerance for resolving the LRO point
Pmin_tol = 1e-9;  % [s] tolerance for resolving the P min point
dT_sup_tol = 1e-12;% [s] tolerance for resolving the point when superheat = 0

N_filter_pts = 10;  % number of points used to low pass 
% filter and get derivatives via backwards differencing

% parameters that may depend on platform (PC/Mac/Linux or laptop/server)
switch computer
    case {'GLNXA64'}
        % on server
        
        addpath('~/refprop/refprop');   % need to add the refprop folders to the path
        t_save = 60*60; % save interval in seconds
        plot_stuff = 1; % whether to plot results when finished
        save_stuff = 1; % whether to save results when finished
        save_periodically = 0;  % whether to save results periodically
        save_parameters_only = 0;   % whether to save everything, or just some key parameters
        plot_periodically = 0;  % whether to plot results periodically
        time_out = 1;   % whether to stop program if it's taking too long
        max_comp_time = 7*60*60;    % [s] how long to wait before timing out
        results_save_dir = pwd; % where to save results to 
           
        
    case {'MACI64','PCWIN64'}
        % on laptop
        
        t_save = 5*60; % save interval in seconds
        plot_stuff = 1;
        save_stuff = 1;
        save_periodically = 1;
        save_parameters_only = 0;
        plot_periodically = 0;
        t_plot = 10; % [s] interval to plot results
        time_out = 1;
        max_comp_time = 6*60*60;
        if strcmp(computer,'PCWIN64')
            results_save_dir = 'C:\Users\jonah\compiled research\model_results';
        else
            results_save_dir = '/Users/jez/School/stanford/compiled research/tank modeling/model_results_data_files';
        end
        
end

%--------------------------------------------------------------------------
% begin initializing things
%--------------------------------------------------------------------------
tic % time the whole code (there's a toc at the end)

% turn off some warnings that are not terribly helpful
% warning('off','MATLAB:nearlySingularMatrix')
warning('off','MATLAB:Axes:NegativeDataInLogAxis')
warning('off','MATLAB:interp2:NaNstrip')

% create a pause button
% but only if we're not on a linux box
if ~strcmp(computer,'GLNXA64')
    button_handle = createButton;
    pause(1e-3);
end

if nargin > 0
    % this script is being called by something else
    % and these parameters are passed
    inputs = varargin{1};
    N_nodes = inputs.N_nodes;
    N_mom = inputs.N_mom;
    rel_tol = inputs.rel_tol;
    constants.C_qdot_lw = inputs.C_qdot_lw;
    constants.C_coalescence = inputs.C_coalescence;
    constants.C_nuc_rate = inputs.C_nuc_rate;
    N_rw = inputs.N_rw;
    
    % turn off some things
    plot_stuff = 0;
    save_stuff = 1;
    save_periodically = 1;
    save_parameters_only = 0;
    plot_periodically = 0;
    time_out = 0;
end


N_ab = N_mom/2; % number of abscissas or weights = # of moments/2

% initialize clocks for saving and plotting
clock_save = clock; 
clock_start = clock;
clock_plot = clock;

% divide the whole tank into nodes or cells
% L_node = L_tank*(1 + 20*alpha_ic)*fill_level/(N_nodes - 0.5);
L_node = L_tank/(N_nodes); % [m] cell/node length
% use N -1 if there are points on the boundaries
V_node = pi*0.25*D_tank^2*L_node;   % [m^3] node/cell volume

% set all nodes to be the same
L_node = L_node*ones(N_nodes,1);
V_node = V_node*ones(N_nodes,1);

% store some things into the constants structure
p = constants.ADQMOM_p;
constants.N_nodes = N_nodes;
constants.N_rw = N_rw;
constants.C_rdot = 1;   % multiplies the rdot expression
constants.C_u_rise = 1; % multiplies the bubble rise rate expression
constants.n_nuc_freq = 3;% exponent on superheat for nucleation rate
constants.C_r_nuc = 1;  % multiplies the size of nucleating bubbles
constants.C_dTs = 1;    % multiplies the superheat
constants.C_x_inj = 1;  % sets the fluid quality at the orifice
constants.V_tank = V_tank;
constants.fluid = fluid;
constants.L_node = L_node;
constants.V_node = V_node;
constants.N_ab = N_ab;

max_L = max_L*L_tank; % rescale from normalized length to actual length [m]

% 6 variables for liquid and vapor
% then N_mom for each node (N_mom/2 abscissas, N_mom/2 weights)
N_dim = 4 + 2*N_rw + N_mom*N_nodes;

% load some fits

% load fit for hesson's data for the orifice mass flux
hesson_fit = load('hesson_fit');
constants.hesson_fit = hesson_fit;

% load PDT table for appropriate fluid
if strcmp(fluid,'N2O')
    PDT = load('N2O_PDT_table','PDT');
    PDT = PDT.PDT;
elseif strcmp(fluid,'CO2')
    PDT = load('CO2_PDT_table','PDT');
    PDT = PDT.PDT;
else
    error('fluid string incorrect. try N2O or CO2')
end

% need this code for using qinterp2 (a fast 2D interpolation script):
Tvec_table = PDT.T;
Pvec_table = PDT.P;

[Tgrid_table, Pgrid_table] = meshgrid(Tvec_table, Pvec_table);
PDT.T = Tgrid_table;
PDT.P = Pgrid_table;

% get ODE solver parameters/constants
switch ode_solver
    case 'ROCK2'
        % orthogonal runge-kutta-chebyshev
        [sigma_R2, gs_term_R2, a_R2, b_R2, c_R2] ...
            = rock2_solver_parameters(ROCK2_stages);
        
        [adaptive, a, b, c, bs, s, p_tilde] = ...
            butcher_tableau(ode_solver, ROCK2_stages, a_R2, b_R2, c_R2);
        % s = ROCK2_stages - 2;
        
        s = ROCK2_stages - 1; % change it now to get the routine right
        
    case 'ROS3P'
        % rosenbrock
        [adaptive, a, b, c, bs, s, p_tilde] = butcher_tableau(ode_solver);
        gamma_ROS = 7.886751345948129e-01;
        
        
    case 'TRBDF2'
        % trapezoidal rule/2nd order backwards difference
        [adaptive, a, b, c, bs, s, p_tilde] = butcher_tableau(ode_solver);
        gamma_TB = a;
        
    otherwise
        % basically a simple runge-kutta type algorithm
        [adaptive, a, b, c, bs, s, p_tilde] = butcher_tableau(ode_solver);
        
end

% set the "adaptive" switch to 0 if we're using constant time step
if force_constant_h
    adaptive = 0;
end

% initialize the error detection flag
constants.error_detected = 0;

% k is a vector used in the ODE integration scheme, initialize it to 0
k_ode = zeros(N_dim,s);

% initialize time
t = 0;

% if not using ADQMOM, set p to 1
if strcmp(ADQMOM, 'off')
    p = 1;
    constants.ADQMOM_p = p;
end

% set the index of the moments that is r^3 (related to volume)
% should = 4 (0,1,2,3) unless using ADQMOM
if p == 1
    constants.V_moment_index = 4;
elseif rem(p*3,1) == 0
    constants.V_moment_index = 1 + 3*p;
end

% store it in the constants structure
V_moment_index = constants.V_moment_index;

% generate initial condition for weights and abscissas
% 2 options here:
% 1 - use gamma distribution to generate moments, then use PD method
% or PD method alternative to turn them into weights and abscissas
% 2 - directly specify the weights and abscissas

% IC_moments = gamma_dist_moments( 5e-10, 1e-3, N_mom, p);
%
% alpha_ic = 1e-12;

% alpha_mom = 4/3 * pi * IC_moments(V_moment_index);
%
% IC_moments = IC_moments * alpha_ic/alpha_mom;
%
% % [r_ic, w_ic] = PD_method(IC_moments);
%
% [r_ic, w_ic] = PD_method_alternative(IC_moments,p);

% based on what r_nuc is initially I would like to make r_ic smaller, but
% code seems to crash if I do
r_ic = logspace(-5, -3, N_ab);
w_ic = ones(1,N_ab);

if ~isreal(r_ic)
    error('imaginary IC')
end

% copy for all nodes
r_ic = repmat(r_ic(:), N_nodes, 1);
w_ic = repmat(w_ic(:), N_nodes, 1);

DQMOM_IC = [w_ic; r_ic.*w_ic];


% reshape so each ROW is one node
w_ic = reshape(w_ic, N_ab, N_nodes)';

% reshape so each ROW is one node
r_ic = reshape(r_ic, N_ab, N_nodes)';

for i = 1:2*N_ab
    mom{1}(:, i) = sum( r_ic.^((i-1)/p) .* w_ic, 2 );
end

% net volume of bubbles, per unit volume of liquid

V_bubi{1} = 4/3*pi*mom{1}(:,V_moment_index);

T_s = Ti; % initialize saturation temp

V_l = fill_level*V_tank; % initialize liquid volume

% get the individual node fill levels
node_level = get_node_levels(V_l, V_bubi{1}, constants.V_node);
guesses.node_level = node_level;

% number of nodes that are full
N_full = sum(node_level == 1);

% total volume of bubbles 
V_bub = sum_over_nodes(V_bubi{1}, node_level, constants.V_node);

% combined volume of bubbles and liquid
V_l_star = V_bub + V_l;

% volume of the ullage (tg = tank gas)
V_tg = V_tank - V_l_star;

% a guess for liquid density is initialized
% it will be needed later when solving nonlinear equation for it
guesses.rho_l = refpropm('D','T',Ti,'Q',0,fluid);

% find the initial pressure
% non-linear equation, have to converge to the actual IC of the system
Pi = converge_on_IC(Ti, V_tg, V_l, V_bubi{1}, PDT, guesses, constants);


% calculate liquid and ullage densities, internal energy of ullage too
% if interpolating a table for properties ('PDT'), use that for the liquid
% density and then use a fit for the saturated vapor properties
% else use refprop and solve nonlinear equations
if strcmp(constants.property_source,'PDT')
    
    rho_l = qinterp2(PDT.T, PDT.P, PDT.D_liq, Ti, Pi/1e3);
    
    % if rho_l is NaN, it means we went outside the bounds of PDT, so
    % instead extrapolate it using interp2 (slower than qinterp2)
    if isnan(rho_l)
        rho_l = interp2(PDT.T, PDT.P, PDT.D_liq, Ti, Pi/1e3, 'spline');
    end
    
    [~, rho_tg, ~, u_tg] = fits_for_getting_P(Pi, fluid);
    
    
elseif strcmp(constants.property_source,'refprop')
    
    rho_l = get_D_from_TP(Ti, Pi, guesses, constants, fluid);
    
    [rho_tg, u_tg] = refpropm('DU','P',Pi/1e3,'Q',1,fluid);
    
end

% initialize pressure
P = Pi;

% initialize saturated vapor density (not the mixture)
rho_tg_sat = rho_tg;

% initialize ullage mass
m_tg = rho_tg * V_tg;

% initialize liquid mass
m_l = rho_l * V_l;

% initialized ullage total internal energy (not specific)
U_tg = m_tg * u_tg;

% initialize liquid temperature
T_l = Ti;

% initialize ullage quality
x_tg = 1;

% initialize wall temperatures
T_lw = Ti*ones(N_rw,1);
T_gw = T_lw;

% initial vector of variabels without the DQMOM variables
y_i_basic = [m_tg; U_tg; T_gw;  m_l; Ti; T_lw];

y_i_DQMOM = DQMOM_IC(:);

% add the DQMOM variables
y = [y_i_basic; y_i_DQMOM];

% load universal constants
[K_b, N_A, h_planck] = universal_constants('boltzmann', 'avagadro', 'planck');

% calculate critical temperature and pressure
[P_cr, T_cr] = refpropm('PT', 'C', 0, '', 0, fluid);
P_cr = P_cr*1e3;

% store things in the constants or guesses structures
constants.D_tank = D_tank;
constants.t_w = t_w;
constants.rho_w = rho_w;
constants.cv_w = cv_w;
constants.Cd = Cd;
constants.A_inj = A_inj;
constants.Po = Po;
constants.T_air = T_air;
constants.h = 0;
constants.k_w = k_w;
constants.K_b = K_b;
constants.P_cr = P_cr;
constants.T_cr = T_cr;
constants.N_A = N_A;
constants.h_planck = h_planck;
constants.g = 9.81;
constants.C_hamaker = 1e-20;
guesses.P = Pi;
guesses.rho_tg = rho_tg;
guesses.rho_l = rho_l;
guesses.Vdot_l = 0;
constants.r_nuc = 0;
guesses.rhodot_tg_sat = 0;

V_bub = 0;
gas_holdup = 0;
deltaT_sup = 0;
Pdot = 0;
rhodot_l = 0;
rhodot_tg = 0;
net_nuc = zeros(N_nodes,1);
f = zeros(N_dim,1);

% flag used to determine if we've already hit the pressure minimum
constants.min_flag = 0;

% the maximum superheat reached so far
constants.deltaT_sup_max = 0;


% parameters for the low pass filter used to calculate some derivatives 
% based on previous time steps using backwards differences

f_cutoff_norm = 1e-1; % supposed to be filter cutoff f / (1/2 sample f)
% so it's basically = 2*dt*f_cutoff
filter_order = 3; % order of the filter

% make a filter with the right order, and set the 3dB down point
filter_handle = fdesign.lowpass('N,F3dB',filter_order, f_cutoff_norm);

% make it a butterworth filter
filter_design = design(filter_handle,'butter');


% turn the node level info into a cell
% this is so that the number of nodes can change
% (only useful when doing adaptive mesh refinement in space)
node_level_cell = cell(1);
node_level_cell{1} = node_level;
node_level = node_level_cell;

h_min_error_count = 0; % counter for the number of times we've had errors 
%from hitting the minimum allowable step size

diff_eqns_error_flag = 0; % error flag for an error in the differential equations
n_increase = 0; % the time step when we last increased the step size
reject_counter = 0; % number of times the current step has been rejected
n_jac = 0; % the time step when we last computed the jacobian
first_save = 1; % flag for whether we're still on the first save 
% (data can be saved periodically, this flag is necessary so the filenames don't get screwed up)
n_coarsened = 0; % the time step we last coarsened the spatial mesh
y_current = y; % current value of the vector of variables
rejected_step = 0; % flag that specifies whether the current step is rejected
running = 1;        % [] switch, 1 = program running, 0 = program stopped
n = 1;              % [] counter

%--------------------------------------------------------------------------
% the main loop -> step through time and solve ODEs
%--------------------------------------------------------------------------

% begin looping
while running == 1;
    
    
    % get dPdt from backwards differences (only used to check if we've
    % reached the P min). bdiff is sort of a legacy piece of code, but it
    % works so leave it be.
    starti = max([n-3, 1]);
    
    Pdot = 0.5*Pdot + 0.5*bdiff(P,starti,n,t,adaptive);
    %     rhodot_l = 0.1*rhodot_l + 0.9*bdiff(rho_l,starti,n,t,adaptive);
    %     rhodot_tg = 0.1*rhodot_tg + 0.9*bdiff(rho_tg,starti,n,t,adaptive);
%     Vdot_l(n+1) = V_tank*bdiff(V_l/V_tank,starti,n,t,adaptive);
    %     Vdot_tg(n+1) = V_tank*bdiff(V_tg/V_tank,starti,n,t,adaptive);
    
    
    
    % filter quantities that are needed from previous time steps to help
    % suppress instability: liquid level and saturated vapor density
    % these will both be differentiated numerically in time
    % (ideally these would be calculated based on current values rather
    % than using previous time steps, but that would turn my whole diffeqns
    % function into a nonlinear equation that would have to be solved)
    
    % only do it if we have enough time steps to get a good calculation
    if n > N_filter_pts && t(n) > 1e-3
        
        % time span of points to filter
        dt_filt = 0.10;
        
        % find the start index
        [~,n_filt_start] = min(abs(t - (t(n) - dt_filt)));
        
        % if there aren't enough points, add more
        if n_filt_start >= n- (filter_order*3 );
            n_filt_start = n - (filter_order*3 + 1);
        end
        
        % compute new number of points to be filtered
        N_filter_pts = n - n_filt_start;
        
        % create a vector of linearly spaced time points that will be used
        % for the filter
        t_filt = linspace(t(n-N_filter_pts),t(n),N_filter_pts+1);
        
        % map the liquid level and saturated vapor density onto the linear
        % time vector
        LL_for_filter = interp1(t(n-N_filter_pts:n),V_l_star(n-N_filter_pts:n)*(L_tank/V_tank), t_filt,'linear');
        rho_tg_sat_for_filter = interp1(t(n-N_filter_pts:n),rho_tg_sat(n-N_filter_pts:n), t_filt,'linear');
        
        % sample frequency
        f_sample = 1/mean(diff(t_filt));
        
        % get the filter cutoff frequency from the sample freq. and the
        % normalized cutoff freq.
        f_cutoff = f_cutoff_norm*f_sample/2;
        
        % only change the filter if the corner frequency got too high or
        % too low (b/c of time step getting small/large)
        if f_cutoff > 100
            % corner freq. too high
            
            % new cutoff frequency
            f_cutoff_norm = 100/(0.5*f_sample);
            
            % if the new cutoff freq. is acceptable, remake the filter
            if f_cutoff_norm < 1
                
                %                 [b_filter,a_filter] = butter(filter_order, f_cutoff_norm, 'low');
                filter_handle = fdesign.lowpass('N,F3dB',filter_order, f_cutoff_norm);
                filter_design = design(filter_handle,'butter');
                
            end
        elseif f_cutoff < 20
            % corner freq. too low
            
            % new cutoff frequency
            f_cutoff_norm = 20/(0.5*f_sample);
            
            % if the new cutoff freq. is acceptable, remake the filter
            if f_cutoff_norm < 1
                %                 [b_filter,a_filter] = butter(filter_order, f_cutoff_norm, 'low');
                filter_handle = fdesign.lowpass('N,F3dB',filter_order, f_cutoff_norm);
                filter_design = design(filter_handle,'butter');
                
            end
        end
        
        % filter the data and calculate the derivatives
        
        %         rho_tg_sat_filtered = filtfilt(b_filter, a_filter, rho_tg_sat);
        rho_tg_sat_filtered = filtfilt(filter_design.sosMatrix, ...
            filter_design.ScaleValues, rho_tg_sat_for_filter);
        
%         guesses.rhodot_tg_sat = bdiff(rho_tg_sat_filtered, 8, 11, t_filt, 0);
        
        guesses.rhodot_tg_sat = finite_diff(rho_tg_sat_filtered, 1, ...
            N_filter_pts+1, N_filter_pts+1, mean(diff(t_filt)), 'backwards', 2);
        
        %         LL_filtered = filtfilt(b_filter, a_filter, LL_for_filter);
        LL_filtered = filtfilt(filter_design.sosMatrix, ...
            filter_design.ScaleValues, LL_for_filter);
        
%         guesses.dLL_dt = bdiff(LL_filtered, 8, 11, t_filt, 0);
        
        guesses.dLL_dt = finite_diff(LL_filtered, 1, ...
            N_filter_pts+1, N_filter_pts+1, mean(diff(t_filt)), 'backwards', 2);
        
        
    else
        % if we couldn't calculate the derivatives, just set to 0
        guesses.rhodot_tg_sat = 0;
        guesses.dLL_dt = 0;
    end
    
    % store the new derivative values
    dLL_dt(n+1) = guesses.dLL_dt;
    rhodot_tg_sat(n+1) = guesses.rhodot_tg_sat;
    
%     Vdot_bub(n+1) = bdiff(V_bub, starti, n, t, adaptive);
    
%     if h > 5*h_min
%         
%         guesses.Vdot_l = 0.5*Vdot_l(n+1) + 0.5*guesses.Vdot_l;
%         
%     else
%         
%         guesses.Vdot_l = 0;
%     end
    
    % store the superheat level in the outer loop (not within the
    % differential equations file) for error checking
    constants.outerloop_superheat = deltaT_sup(n);
    
    % check for P min
    if (t(n) > 0.1) && (Pdot > 10)
        if constants.min_flag == 0;
            constants.min_flag = 1;
            disp('P min')
            t_min = t(n);
            n_min = n;
        end
    end
    
    % if the current superheat is higher than the recorded max, update
    if deltaT_sup(n) > constants.deltaT_sup_max
        constants.deltaT_sup_max = deltaT_sup(n);
    end
    
    % current value of time and time step (from last succesful step)
    dt = t(n) - t(max([1, n-1]));
    constants.t = t(n);
    constants.dt = dt;
    
    
    % check if the pressure got to low, stop if it did
    if P(end) < 2e6
        running = 0;
        stop_reason = 'P got below 2 MPa';
    end
    
    
    
    % if I'm just playing around, print status at each step
    if nargin == 0
        
        fprintf(['t = %#4.4g, dt = %#4.4g, P = %#4.4g, alpha = %#4.4g, dT_sup = %#6.6g,' ...
            ' V_bub = %#4.4g, T_l = %#4.4g, m_l/m_li = %#4.4g,'...
            ' m_tg/m_tgi -1 = %#4.4g, fill_level%% = %#4.4g\n'],...
            t(n), t(n) - t(max([1, n-1])), P(n)/1e6, gas_holdup(n), deltaT_sup(n), ...
            V_bub(n), T_l(n), m_l(n)/m_l(1), m_tg(n)/m_tg(1)-1, 100*fill_level(n));
        
        
    end
    
    % if we're using an adaptive scheme, check for some things and refine
    % the time step accordingly
    if  n > 1 && adaptive == 1
        
        % --- check if we're close to LRO ---
        
        % slope of fill level curve
        LRO_slope = bdiff(fill_level,starti,n,t,adaptive);
        
        % projected t_LRO
        t_LRO = -fill_level(n)/LRO_slope + t(n);
        
        h_LRO = t_LRO - t(n); % distance to t_LRO
        
        % if the step we're about to take is >3/4 the distance to LRO
        % and the distance te LRO is bigger than the tolerance
        if (h > 2*h_LRO && h_LRO > LRO_tol) && (h_LRO > 0);
            
            % set h to 1/2 the distance to LRO (ie refine)
            h = 0.5*h_LRO;
            
        end
        
        % refine as we get close (LRO_slope < 0)
        if (LRO_slope*h < -0.003/100) && (fill_level(n) < 0.1/100)
            h = h/4;
        elseif (LRO_slope*h < -0.03/100) && (fill_level(n) < 1/100);
            h = h/2;
        elseif (LRO_slope*h < -0.3/100) && (fill_level(n) < 5/100);
            h = h/2;
        end
        
        
        
        % also check if we're close to going subcooled -> superheated
        % really only was useful when I wasn't handling the initial
        % conditions well and it started subcooled. Now it's not really
        % useful
        
        % slope of dT_superheat curve
        dTs_slope = bdiff(deltaT_sup,starti,n,t,adaptive);
        
        % projected t_sup (t when dT_sup = 0)
        t_sup = -deltaT_sup(n)/dTs_slope + t(n);
        
        
        h_sup = t_sup - t(n); % distance to t_sup
        
        % if the step we're about to take is >3/4 the distance to the
        % crossover point and the distance is bigger than the tolerance
        if (dTs_slope > 0 && deltaT_sup(n) < 0) || (dTs_slope < 0 && deltaT_sup(n) > 0 )
            
            if (h > 0.5*h_sup && h_sup > dT_sup_tol) && (h_sup > 0);
                
                % set h to 1/2 the distance to crossover point (ie refine)
                h = 0.5*h_sup;
                disp('refining based on superheat being near 0')
                
            end
            
        end
        
        % hysteresis in the bubble nucleation process
        % if we're using it, we check how close we are to Pmin and refine
        % note: doesn't seem to work that well
        
        if constants.include_hysteresis
        % also check if we're close to Pmin (it's been crashing there...)
            if constants.min_flag == 0 && t(n) > 0.15
                % slope of Pdot curve
                k = 1;
                clear Pdot2
                for j = n-25:n % use last 25 points
                    Pdot2(k) = bdiff(P,j-3,j,t,adaptive); % use backward difference, 3rd order
                    k = k + 1;
                end

                % fit line
                fit_line_data = [ones(k-1,1) t(n-(k-2):n)']\Pdot2';
                
                % find point when Pdot = 0
                t_Pd = -fit_line_data(1)/fit_line_data(2);

                h_Pd = t_Pd - t(n); % distance to Pdot = 0 point

                if (h > 0.01 * h_Pd && h > Pmin_tol) && h_Pd > 0
                    h = h_Pd/5;
                    disp('refining because close to Pmin')
                end

            end
        end
        
        
    end
    
    
    
    if strcmp(ode_solver,'ROS3P')
        % if using rosenbrock scheme, calculate jacobian
        % I might have to move this inside error_OK while loop if I do adaptive
        % mesh refinement...
        
        %         fy = f(t(i), y(:,i));
        
        % set the constant for step number to 1
        constants.step = 1;
        
        [fy, debug_data] = diffeqns(y_current, constants, guesses, PDT); % get the nominal value of the diff eqns
        
        dfdy = zeros(N_dim,N_dim); % initialize the jacobian to zeros
        
		% loop through and calculate partial derivative for each dimension
		% computes one sided derivative (ie (y(x+dx) - y(x))/dx)
        for k = 1:N_dim
            dy = 1e-5; % relative change in the variable
            dy = dy*abs(y_current(k)); % absolute change in the variable (make sure it's positive)
            if dy == 0 % check if the derivative is zero, if so set it to a small number
                dy = 1e-9;
            end
            y_plus = y_current(:); % initialize 
            y_plus(k) = y_current(k) + dy; % increment the kth value
            
            [fy_plus, debug_data] = diffeqns(y_plus, constants, guesses, PDT); % calculate new value of derivatives
            
            dfdy(:,k) = ( fy_plus - fy )/dy; % insert column of jacobian
        end
        
        
        
    end
    
	% while loop that keeps going until we've satisfied error requirements
	% solves differential equations for a time step
    error_OK = 0;
    while error_OK == 0
        
			% insert this pause only if on linux
			% purpose: to allow pause button to work
            if ~strcmp(computer,'GLNXA64')
                pause(1e-6)
             end  
			 
        % solving differential equations
        
		% total number of dimensions
		% (this changes if I do adaptive mesh refinement)
        N_dim = 4 + 2*N_rw + N_mom*N_nodes;
        
        constants.h = h; % set the constant to include current time step size
        
        error_flag = 0; % initialize error flag
        
        % clear this variable (if it exists)
		% (important for adaptive mesh refinement -> size of variable changes)
        if exist('k_ode','var')
            clear k_ode
        end
        
		% depending on the scheme we're using (implicit vs explicit), do different things
        if strcmp(ode_solver,'TRBDF2')
            % implicit method!
            
            constants.step = 1; % set the constant so the diff eqns think we're on the first step (of a time step)
            
            
			% call diffeqns and compute current value
			% use a try/catch here to deal with errors within diffeqns
            try
                [fy, debug_data] = diffeqns(y_current, constants, guesses, PDT);
                if debug_data.diff_eqns_error_flag % check for error flag within diffeqns
                    error_flag = 1;
                    disp('error_flag tripped in diffeqns')
                    
                end
            catch
                disp('threw an error calling diffeqns, beginning of TRBDF2')
                error_flag = 1;
                
            end
            
            % if necessary, compute jacobian
			% it's computationally intensive so only do it if needed
            
			% conditions that suggest we need a new jacobian
			% the time step is important because we're actually computing a quantity that involves the jacobian and the time step, not just the jacobian
            new_jacobian_conditions = [ n == 1; % it's the first time step
                (n>1 && rejected_step) && ( n~=n_jac && reject_counter>2); % it's not the first step, we rejected the last step, we didn't already just compute it, and we've rejected the step more than 2x
                (n > 1) && ( (abs(h - h_jac)/h_jac > 10) || (h/h_jac < 0.1) ); % it's not the first step, the time step has grown or shrunk by 10x since we computed it last
                (n > 1) && ( (n - n_jac) > 100 )]; % it's not the first step and it's been 100 steps since we last calculated it
            
			% if we've violated one of the above, compute new jacobian
            if sum(new_jacobian_conditions) > 0
                disp('computing new jacobian')
                
				% store the time step and step size
                h_jac = h;
                n_jac = n;
                
				% if fac exists, clear it
				% fac is used by Matlab's numjac routine to store some information about the jacobian
                if ~exist('fac','var')
                    fac = [];
                end
                
                    clear i_w_empty i_g_empty % clear variables
					
					% find the indices for weights and abscissas associated with empty nodes%
					% don't want to bother computing jacobian at these points!
                    i_w_empty = (4+2*N_rw) + [(N_full(n)+2)*N_ab:N_nodes*N_ab];
                    i_g_empty = (4+2*N_rw) + N_nodes*N_ab + [(N_full(n)+2)*N_ab:N_nodes*N_ab];
                    i_empty = [i_w_empty(:); i_g_empty(:)];
                
					% if we're using matlab's numjac routine, use that to calculate df/dy
                    if constants.use_numjac
						% I had to modify it slightly (hence it's now numjac_no_t.m) to deal with the fact that it was expecting a time variable as well as y (f(t,y)).
						[dfdy, fac] = numjac_no_t(@(y) diffeqns(y, constants, guesses, PDT), y_current, fy, [], abs(y_current)*1e-8,fac,0);
					
					% if we're not using numjac, just loop through and do a central difference
					else
						
						dfdy = zeros(N_dim,N_dim); % initialize
                
						% calculate each column (1 for each dimension of y)
						for k = 1:N_dim
				
							% check if current dimension is one of an empty node
							% if so, just set it to zeros
							if sum(i_empty == k) > 0
								% it's empty
								dfdy(:,k) = zeros(N_dim,1);
							else
								
								dy = 1e-4; % set relative change in y
								dy = dy*abs(y_current(k)); % absolute change in y
								if dy == 0
									dy = 1e-9; % if y is 0, just use a small number
								end
								y_plus = y_current(:); % initialize y + dy
								y_plus(k) = y_current(k) + dy; % add in the incremented variable
								constants.step = 0; % set this to 0 so diffeqns knows what's going on
								[fy_plus] = diffeqns(y_plus, constants, guesses, PDT); % calculate f at y + dy
								
								% repeat for y - dy
								y_minus = y_current(:);
								y_minus(k) = y_current(k) - dy;
								constants.step = 0;
								[fy_minus] = diffeqns(y_minus, constants, guesses, PDT);
								
								
								% compute central difference and add the column to the jacobian
								dfdy(:,k) = (fy_plus - fy_minus)/(2*dy);
							end
						end
                    end
					
				% identity matrix the right size
                I = eye(N_dim);
                
				% see equation 5.89 in my thesis
				% this is the jacobian for the numerical method
				% (gamma_TB defines the relative size of the trapezoidal rule step and the 2nd order bacwards difference step, = 2-sqrt(2))
                J = I - gamma_TB*h/2 * dfdy;
                
				% I played around with different methods for computing the inverse of J because it's often ill-conditioned
				% but ended up sticking with inv(J) because it worked well enough
%                 if rcond(J) < 1e-16
%                     % poorly conditioned!
%                     inv_J = pinv(J);
%                 else
                    inv_J = inv(J);
%                 end
            end
            
            
            % initialize x, F, d for TR
			% I use somewhat different notation here because that's what was in the papers I was using for info on TRBDF2
            % u_n is u^n in my thesis and is y at the current step
			% x_k is a proxy for u -> it starts as equal to u_n, but then as I iterate to find u at n+gamma, I store the intermediate value as x (k is the counter for the iteration in Broyden's method)
			u_n = y_current;
            x_k = y_current;
            F_k = -gamma_TB * h * fy; % F_k = F(u_k)
            
            d_k = -inv_J * F_k; % delta x needed from k to k+1
            x_kp1 = x_k + d_k; % x at k+1
            k = 1; % initialize k
            
            % trapezoidal rule step
			% keep iterating until converged or if we determine it's not converging
            not_converged = 1;
            while not_converged
                
                % evaluate function at current point x_kp1
                try
                    constants.step = 2;
                    [f_kp1, debug_data] = diffeqns(x_kp1, constants, guesses, PDT);
                    if debug_data.diff_eqns_error_flag
                        error_flag = 1;
                        disp('error_flag tripped in diffeqns')
                        break
                    end
                catch
                    disp('threw an error calling diffeqns, TR')
                    error_flag = 1;
                    break
                end
                
				% F at k + 1
                F_kp1 = x_kp1 - y_current - gamma_TB*h/2*( f_kp1 + fy );
                
				% new value of u
                u_k = inv_J * F_kp1;
                
				% just storage to simplify next expression
                c_k = d_k' * (d_k + u_k);
                
				% new inverse of J (eqn 5.93)
                inv_J = inv_J - 1/c_k * ( u_k * d_k' ) * inv_J;
                
				% increment k
                k = k + 1;
                
				% store old value of k+1 in k
                x_k = x_kp1;
                
                F_k = F_kp1;
                
                d_k_old = d_k;
                
				% new delta x
                d_k = -inv_J * F_k;
                
				% new x at k+1
                x_kp1 = x_k + d_k;
                
				% relative change in r (eqn 5.95)
                r = norm(d_k)/norm(d_k_old);
                
				% if r is small enough, we've converged (5.96)
                if k > 1 && r/(1-r)*norm(d_k./u_n) < newton_tol*rel_tol
                    not_converged = 0; % set flag to exit loop
                end
                
				% if we've had too many itertions, set error flag and exit loop
                if k > max_iter
                    disp('not converging in TR')
                    error_flag = 1;
                    break
                end
                
            end
            
            if error_flag ~= 1
                % don't bother if we already have an error
                % initialize x, F, d for BDF2
                
                u_npg = x_kp1; % value of u at n + gamma
                
				% get f(u at n+gamma)
                try
                    constants.step = 3;
                    [f_npg, debug_data] = diffeqns(x_kp1, constants, guesses, PDT);
                    if debug_data.diff_eqns_error_flag
                        error_flag = 1;
                        disp('error_flag tripped in diffeqns')
                        
                    end
                catch
                    disp('threw an error calling diffeqns, between TR and BDF2')
                    error_flag = 1;
                    
                end
                
				% rename variables
                x_k = u_npg;
                
				% equation (5.90)
                F_k = u_npg - (1 - gamma_TB)/(2-gamma_TB) * h * f_npg ...
                    - 1/(gamma_TB*(2-gamma_TB))*u_npg +...
                    (1-gamma_TB)^2/(gamma_TB*(2-gamma_TB)) * u_n;
                
                d_k = -inv_J * F_k;
                x_kp1 = x_k + d_k;
                k = 1;
                
                % BDF2 (similar to above)
                
                not_converged = 1;
                while not_converged
                    
                    % evaluate function at current point x_kp1
                    try
                        constants.step = 4;
                        [f_kp1, debug_data] = diffeqns(x_kp1, constants, guesses, PDT);
                        if debug_data.diff_eqns_error_flag
                            error_flag = 1;
                            disp('error_flag tripped in diffeqns')
                            break
                        end
                    catch
                        disp('threw an error calling diffeqns, BDF2')
                        error_flag = 1;
                        break
                    end
                    
                    F_kp1 = x_kp1 - (1 - gamma_TB)/(2-gamma_TB) * h * f_kp1 ...
                        - 1/(gamma_TB*(2-gamma_TB))*u_npg +...
                        (1-gamma_TB)^2/(gamma_TB*(2-gamma_TB)) * u_n;
                    
                    u_k = inv_J * F_kp1;
                    
                    c_k = d_k' * (d_k + u_k);
                    
                    inv_J = inv_J - 1/c_k * ( u_k * d_k' ) * inv_J;
                    
                    k = k + 1;
                    
                    x_k = x_kp1;
                    
                    F_k = F_kp1;
                    
                    d_k_old = d_k;
                    
                    d_k = -inv_J * F_k;
                    
                    x_kp1 = x_k + d_k;
                    
                    r = norm(d_k)/norm(d_k_old);
                    
                    if k > 1 && r/(1-r)*norm(d_k./u_npg) < newton_tol*rel_tol
                        not_converged = 0;
                    end
                    
                    if k > max_iter
                        disp('not converging in BDF2')
                        error_flag = 1;
                        break
                    end
                    
                end
                
				% now have values at n + 1
                u_np1 = x_kp1;
                y_new = u_np1;
                
				% compute error estimate (5.9)
                error_estimate = (bs(1) - b(1))*u_n + ...
                    (bs(2) - b(2))*u_npg + ...
                    (bs(3) - b(3))*u_np1;
                
                g_star = y_new - error_estimate; % this is basically the value of y predicted by a lower order integrator
				% it's needed because some of the later code was written to take the two different estimates and compare them to get the error, rather than taking the error directly. So this just maintains compatibility
                
            else
                % we threw an error
                error_flag = 1;
                g_star = y_new;
            end
            
        else
            % NOT TR-BDF2 escheme
            
            if strcmp(ode_solver,'ROS3P')
			% 3rd order rosenbrock scheme
			% need to look in lab notebook for info on this to comment it
                
                I = eye(N_dim); % identity matrix
                
                k_ode = zeros(N_dim,3); % initialize
                
                E_mat = I/(h*gamma_ROS) - dfdy;
                
                k_ode(:,1) = E_mat\fy;
                
                for i = 2:3
                    a_k_term = sum( (ones(N_dim,1) * a(i,1:i-1)) .* k_ode(:,1:i-1) , 2 );
                    c_k_term = sum( (ones(N_dim,1) * c(i,1:i-1)) .* k_ode(:,1:i-1) , 2 );
                    
                    y_intermediate = y_current + a_k_term; % intermediate value of y
                    constants.step = i;
                    [fy, debug_data] = diffeqns(y_intermediate, constants, guesses, PDT);
                    
                    
                    RHS = fy + 1/h * c_k_term;
                    
                    k_ode(:,i) = E_mat\RHS;
                end
            else
                % not ROS3P or TRBDF2 -> we're using a traditional explicit scheme
                
                % what we're basically doing (doesn't have to be 6 stages)
                %         k1 = h*diffeqns(y(:,n));
                %         k2 = h*diffeqns(y(:,n) + a(2,1)*k1);
                %         k3 = h*diffeqns(y(:,n) + a(3,1)*k1 + a(3,2)*k2);
                %         k4 = h*diffeqns(y(:,n) + a(4,1)*k1 + a(4,2)*k2 + a(4,3)*k3);
                %         k5 = h*diffeqns(y(:,n) + a(5,1)*k1 + a(5,2)*k2 + a(5,3)*k3 + a(5,4)*k4);
                %         k6 = h*diffeqns(y(:,n) + a(6,1)*k1 + a(6,2)*k2 + a(6,3)*k3 + a(6,4)*k4 + a(6,5)*k5);
                %
                %         k = [k1, k2, k3, k4, k5, k6];
                
                for i = 1:s
                    % s = number of stages in the scheme
                    
                    if i == 1
                        % first stage of the scheme
                        % f = f( t(n) , y(n) )
                        
                        if  (n > 1 && ~isnan(f_np1(1))) && (length(f_np1) == N_dim)
                            % we calculated it at the end of last time step
                            % and we haven't refined in space since then
                            f = f_np1;
                        else
                            % we haven't had a time step yet!
                            constants.step = 1;
                            try
                                [f, debug_data] = diffeqns(y_current, constants, guesses, PDT);
                                if debug_data.diff_eqns_error_flag
                                    error_flag = 1;
                                    disp('error_flag tripped in diffeqns')
                                end
                            catch
                                disp('threw an error calling diffeqns, first step')
                                error_flag = 1;
                            end
                        end
                    else
                        % not the first stage of the scheme
                        % f for k(2) = f( t(n) + c(2)*h , y(n) + a(2,1)*k(1) )
                        % f for k(3) = f( t(n) + c(3)*h , y(n) + a(3,1)*k(1) + a(3,2)*k(2) )
                        % and so on
                        constants.step = i;
                        
                        a_k_term = sum( (ones(N_dim,1)*a(i,1:i-1)).*k_ode(:,1:i-1) ,2 );
                        
                        y_intermediate = y_current + a_k_term;
                        
                        variables = unpack_y(y_intermediate, constants);
                        
                        T_l_new = variables.T_l;
                        g_q_new = variables.g_q;
                        w_q_new = variables.w_q;
                        
                        r_q_new = g_q_new./w_q_new;
                        
                        V_bubi_new = 4/3*pi*sum( r_q_new.^((4-1)/p) .* w_q_new, 2 );
                        
                        % check to see if we're taking a step that'll cause errors
                        error_conditions = [isnan(sum(y_intermediate(:))), (T_l_new < 220), (T_l_new > 305), ...
                            (T_l_new > (T_l(n) + 0.25)),  (T_l_new < (T_l(n) - 5)),...
                            (max(r_q_new(:)) > 0.2), (max(V_bubi_new) > 1)];
                        
						% if there are errors, set error flag and break
                        if sum(error_conditions) > 0
                            disp('we''re taking a bad step')
                            f = ones(size(f));
                            error_flag = 1;
                            % need to fill out k_ode otherwise the dimensions will
                            % be wrong and I'll get another error
                            k_ode = ones(N_dim, s);
                            break
                        else
                            % no errors
							
							% check for negative values in y
                            [ind_neg, min_val] = min(y_intermediate);
                            if min_val < 0
                                fprintf('negative part of y. index: %0.d\n', ind_neg)
                            end
                            
							% compute derivatives at y_intermediate
                            try
                                [f, debug_data] = diffeqns(y_intermediate, ...
                                    constants, guesses, PDT);
                                if debug_data.diff_eqns_error_flag
                                    error_flag = 1;
                                    disp(['error_flag tripped in diffeqns, i = ' num2str(i)])
                                    
                                end
                            catch ME
								% if there was an error in diffeqns, report it to the screen and set error flag
                                disp(['threw an error calling diffeqns, i = ' num2str(i)])
                                disp( getReport(ME))
                                error_flag = 1;
                            end
                            
                        end
                    end
                    
                    % the new value of k (f is dy/dt)
                    k_ode(:,i) = f*h;
                    
                end
                
            end
            
            if strcmp(ode_solver, 'ROCK2')
                % ROCK2 solver - computing y_new is a little different
                
                % already calculated g_s-2
                g_sm2 = y_intermediate;
                
                % s = s - 1
                
                f_g_sm2 = f;
                g_sm1 = g_sm2 + h * sigma_R2 * f_g_sm2;
                
                % s = s
                
                
                try
                    [f_g_sm1, debug_data] = diffeqns(g_sm1, ...
                        constants, guesses, PDT);
                    if debug_data.diff_eqns_error_flag
                        error_flag = 1;
                        disp(['error_flag tripped in diffeqns, 1st finishing stage'])
                        
                    end
                catch ME
                    disp('threw an error calling diffeqns, first of the finishing stages')
                    disp( getReport(ME))
                    error_flag = 1;
                end
                
                %         f_g_sm1 = f(t(i),g_sm1);
                g_star = g_sm1 + h * sigma_R2 * f_g_sm1;
                
                y_new = g_star + h * gs_term_R2 * (f_g_sm1 - f_g_sm2);
                
                
            else
                % runge-kutta type method
                
				% error checking
                if length(k_ode(1,:)) ~= length(b)
                    disp('dimensions are wrong')
                end
                
				% update y_new
                y_new = y_current + (k_ode*b);
                
            end
            
        end
        
        if error_flag == 0
            % need to check the new y to make sure it's OK too
            
            variables = unpack_y(y_new, constants);
            
            T_l_new = variables.T_l;
            g_q_new = variables.g_q;
            w_q_new = variables.w_q;
            
            r_q_new = g_q_new./w_q_new;
            
            % check to see if we're taking a step that'll cause errors
            error_conditions = [isnan(sum(y_new(:))), (T_l_new < 220), (T_l_new > 305), ...
                (T_l_new > (T_l(n) + 0.25)),  (T_l_new < (T_l(n) - 5)),...
                (max(r_q_new(:)) > 0.2)];%, (max(V_bubi_new) > 1)];
            % (min(y_new) < 0),
            
            if sum(error_conditions) > 0
                disp('we''re taking a bad step')
                f = ones(size(f));
                error_flag = 1;
            end
        end
        
        
        % ----- check error and figure out new step size -----
        
        if adaptive == 1
            % using adaptive scheme, need to check error and pick new time
            % step
            
            if strcmp(ode_solver, 'ROCK2') || strcmp(ode_solver, 'TRBDF2')
                % error calculated a little differently than runge-kutta
                % type solvers
                err = y_new - g_star;
            else
                % runge-kutta type solver
                err = k_ode*(b - bs);   % absolute error (diff. between 5th and 4th order estimates of y(n+1) - y(n))
            end
            
            
            % calculate the relative error in each variable
			% the commented out code is from different methods I've used to calculate this
            for j = 1:N_dim
%                 if exist('i_empty','var') && ~any(j == i_empty)
                    rel_err(j) = abs( err(j) /( mean([y_current(j) y_new(j)])  + 1e-8) );  % relative error
%                 else
%                     rel_err(j) = 1e-16;
%                 end

% if j > 4+2*N_rw
%     if y_new(j) > 1e6
%         rel_err(j) = abs( (log(y_new(j)) - log(g_star(j)))/( log(y_current(j)) + 1e-8) );
%     end
% end
            end
            
            
            [~, ind_max_rel_err(n+1)] = max(rel_err);  % fix rel_err to the maximum finite value of rel_err
            
			% calculate relative error (scalar) based on the desired norm
            switch constants.error_norm
                case 'L-infinity'
                    rel_err = max(rel_err);  % fix rel_err to the maximum finite value of rel_err
                case 'L-2'
                    rel_err = sqrt( sum( rel_err.^2 ) );
                case 'L-1'
                    rel_err = sum( rel_err );
                otherwise
                    error('incorrect constants.error_norm expression. try ''L-infinity'' or L-2 or L-1')
            end
            
            
            abs_err = abs(err); % absolute error
            
            abs_err = max(abs_err(isfinite(abs_err)));  % do the same for abs_err
            
            
            
            % check for possible problems: isempty statements are in case
            % abs and rel err are both full of non-finite values
            % isnan checks for nan's
            % isreal checks for imaginary numbers
            error_conditions2 = isempty(rel_err) + ...
                isempty(abs_err) +  ...
                isnan(sum(err)) + ...
                ~isreal(sum(y_new)) + ...
                (deltaT_sup(n) < 0) + ...
                error_flag;
            
            % if any of those fail, set rel_err large so that the step gets
            % recomuputed
            if error_conditions2 > 0
                rel_err = 1;
                disp('encountered a problem with the error terms')
            end
            
            if isempty(abs_err)
                abs_err = 1;
            end
            
            % pick new step size
            
            error_estimate = rel_err;
            
			% if it's the first step or it's been a long time,
			% reset the error terms for step size control (ie delete history)
            if n == 1 || (n > n_reset + 50)
                h_past = h; % previous step size
                error_estimate_past = error_estimate; % previous step error estimate
                sh_max = constants.sh_max; % max possible factor to increase step size by
                n_reset = n; % time step we last reset
            end
            
			% sh is the factor by which the time step is increased/decreased
			% eg: h (new) = h (old) * sh
			% it's basically picked to keep the error close to the tolerance (but smaller than)
			
			% previous way of calculating it
            % sh = 0.7 * (rel_tol/rel_err).^(1/(p_tilde + 1));
            
			% new way
            sh = (rel_tol/error_estimate)^( 1/p_tilde );
            
			% if it's not the first step and we haven't rejected the step
			% compare to previous value of sh and pick minimum of the 2
            if n > 1 && rejected_step ~= 1
                
                sh_past = (error_estimate_past/rel_tol)^(1/p_tilde) * sh^p_tilde * (h / h_past);
                
                sh = min(sh, sh_past);
                
            end
            
			% this makes sure sh isn't too small or too big
			% and also multiplies it by 0.8
            sh = min( sh_max, max( 0.1, 0.8 * sh) );
			
			
			% if error is acceptable or if the step size is below the min,
			% accept the step
            if ( rel_err < rel_tol && abs_err < abs_tol) || (h < h_min)

				% error is OK!
                error_OK = 1;
                
				% we didn't reject the step
                rejected_step = 0;
                reject_counter = 0;
                
                % reset the max allowable sh
                sh_max = constants.sh_max;
				
xxxxx				% if the last step was rejected, keep sh from being > 1
				% this doesn't make sense because rejected_step is set to 0 a few lines above
                if rejected_step == 1
                    % we rejected the last step
                    sh = min(sh, 1);
                end
                
				% store h as past h
                h_past = h;
				
				% store error estimate as past error estimate
                error_estimate_past = error_estimate;
                
                % if it's not the first step and we're close to running out of liquid,
				% stop the train because we want to get off now (we're done!)
                if ((n > 1) && ((h_LRO < LRO_tol) && (h_LRO > 0))) && (fill_level(n) < 0.01)
                    % distance to LRO is less than LRO_tol
                    running = 0;
                end
                
				% if step size is too small,
				% alert the user, reset sh, and increment a counter
				% if we've had too many steps with time step small, stop the program
                if h < h_min
                    fprintf('h got too small. exceeded tolerance by %6.4g%%\n',100*rel_err/rel_tol);
                    h_min_error_count = h_min_error_count + 1;
                    sh = 2*h_min/h;
                    if h_min_error_count > 200
                        running = 0;
                        stop_reason = 'exceeded h_min too many times';
                    end
                end
                
            else
				% current step not meeting error requirements
				
				% we're rejecting the step
                rejected_step = 1;
				
				% increment a counter for how many times we've rejected THIS step
                reject_counter = reject_counter + 1;
                
                % figure out where error is (not really unpacking y, just a poor function name)
                unpack_y(y_new, constants, ind_max_rel_err(n+1), rel_err);
                
				% print error info to screen for debugging
                %                 fprintf(['max rel err = %6.4g, ind of max rel err = %6.4g\n'...
                %                     'err(ind_max_rel_err) = %8.6g, y(ind_max_rel_err,n+1) = '...
                %                     '%8.6g, y(ind_max_rel_err,n) = %8.6g\n'], ...
                %                     rel_err, ind_max_rel_err, err(ind_max_rel_err), ...
                %                     y(ind_max_rel_err, (n+1):-1:n))
                
                
                % if for some reason my error terms are exactly 0,
				% it means something weird happened and we should reduce step size
                if rel_err == 0 || abs_err == 0
                    % something odd happened, so reduce step size a lot
                    disp('something weird happened. reducing step size')
                    sh = sh_min;
                    
                else % otherwise, check for other weird things
                    if rel_err == 1
                        % error flag was probably tripped, or one of the
                        % other error conditions
                        sh = sh_min;
                    end
                    
					% no weird things -> decrease sh appropriately
                    sh = max([0.8 * sh, sh_min]);
					
					% reset sh_max so we don't increase step size too much right after rejecting a step
                    sh_max = 1;
                    
                    disp('rejected step')
                    
                    % commented out code below is from my old method of updating the time step
                    
                    %                     if rel_err/rel_tol > abs_err/abs_tol
                    %                         % if relative error is a bigger problem than
                    %                         % absolute, update step size based on relative
                    %                         % error. Else, use absolute
                    %
                    %                         % refinement algorithm from Numerical
                    %                         % Analysis, by Burden, section 5.5
                    %                         disp('refining step size')
                    %
                    %                         sh = 0.84*( rel_tol*h / (2*rel_err) )^(1/4);
                    %                     else
                    %                         sh = 0.84*( abs_tol*h / (2*abs_err) )^(1/4);
                    %                     end
                    
                end
                
                %                 if sh < 0.1
                %                     % if it looks like the step size would be reduced too
                %                     % much, only reduce it by 1/10
                %                     sh = 0.1;
                %                 elseif sh > 4.0
                %                     % similarly if it's too big, only make it 4x bigger
                %                     sh = 4.0;
                %                 end
                %
                %                 % update step size
                %                 h = h*sh;
                %
                %                 % minimum step size set by computer's precision
                %                 h_min = 16*eps(t(n));
                %
                %                 % self explanatory I think
                %                 if h > h_max
                %                     h = h_max;
                %                 elseif h < h_min
                %                     h = h_min;
                %                 end
                
            end
            
			% update time step
            h = h*sh;
            
			% deal with adaptive mesh refinement if it's turned on
            if constants.adaptive_mesh_refinement
                
                % initialize error control flag
                spatial_error_OK = 1;
                
                if error_OK == 1
                    % only bother to refine in space if we're taking an OK
                    % step already
                    
                    % first check spatial error
                    points_to_refine = 0;
                    points_to_coarsen = 0;
                    
                    % need to be careful about the points at the liquid level
                    % w and g go to zero on the other side. using N_full from
                    % last step
                    
                    % loop through and find relative change between two points
                    % (relative to the interface between them)
                    % then if the diff is too high, tag for refinement
                    for k = 1:N_full(n)
                        
						% change in w and g
                        dw = w_q_new(k+1,:) - w_q_new(k,:);
                        dg = g_q_new(k+1,:) - g_q_new(k,:);
                        
						% mean value of w
                        w_bar = w_q_new(k,:) + (w_q_new(k+1,:) - w_q_new(k,:))/...
                            (L_node(k+1) + L_node(k)) * L_node(k);
                        
						% mean value of g
                        g_bar = g_q_new(k,:) + (g_q_new(k+1,:) - g_q_new(k,:))/...
                            (L_node(k+1) + L_node(k)) * L_node(k);
                        
						% relative change in w and g
                        rel_dw = dw./w_bar;
                        rel_dg = dg./g_bar;
                        
						% max change in relative w and g
                        max_delta_rel = max( abs([ rel_dw(:) rel_dg(:)]) );
                        
						% if the largest difference exceeds the tolerance
						% and we haven't already tagged this point for refinement
						% then tag the point for refinement
                        if max_delta_rel > max_spatial_rel_delta_tol
                            if points_to_refine(end) ~= k
                                points_to_refine = [points_to_refine; k];
                            end
                        end
                        
						% I think this is old code for coarsening the mesh (current code only refines)
                        %                     if n > n_coarsened + 3
                        %                         if max_delta_rel < min_spatial_rel_delta_tol
                        %                             if length(points_to_coarsen) < 2 && L_node(k) < max_L
                        %                                 if points_to_coarsen(end) ~= k
                        %                                     points_to_coarsen = [points_to_coarsen; k];
                        %                                 end
                        %                             end
                        %                         end
                        %                     end
                    end
                    
					% if we tagged points for refinement, create the new point and update appropriate variables
                    if length(points_to_refine) > 1
                        
						% remove the zero
                        points_to_refine = points_to_refine(2:end);
						
						% initialize error flags
                        spatial_error_OK = 0;
                        error_OK = 0;
                        
                        % use old points to refine space, not new points
                        variables = unpack_y(y_current, constants);
                        
                        g_q_new = variables.g_q;
                        w_q_new = variables.w_q;
                        
                        y_old = y_current;
                        L_node_old = L_node;
                        V_node_old = V_node;
                        N_nodes_old = N_nodes;
                        
                        for j = 1:length(points_to_refine)
                            
                            % the point to refine (between k and k+1)
                            % need to correct for points already added -> j-1
                            k = points_to_refine(j) + j - 1;
                            
                            if k > N_full(n) - 1
                                disp('Warning! Trying to refine point near the liquid level.')
                            end
                            
                            % calculate values at new points
                            w_bar(1,:) = w_q_new(k,:) + (w_q_new(k+1,:) - w_q_new(k,:))./...
                                (L_node(k+1) + L_node(k)) * mean(L_node(k:(k+1)));
                            
                            g_bar(1,:) = g_q_new(k,:) + (g_q_new(k+1,:) - g_q_new(k,:))/...
                                (L_node(k+1) + L_node(k)) * mean(L_node(k:(k+1)));
                            
                            if (L_node(k) >= L_node(k+1)) || (k > N_full(n) - 1)
                                % divide k in half
                                
                                % values at the new point
                                w_star = w_q_new(k,:) + (w_q_new(k+1,:) - w_q_new(k,:))./...
                                    (L_node(k+1)/2 + L_node(k)/2) * L_node(k)*0.25;
                                
                                g_star = g_q_new(k,:) + (g_q_new(k+1,:) - g_q_new(k,:))/...
                                    (L_node(k+1)/2 + L_node(k)/2) * L_node(k)*0.25;
                                
                                % values at shifted k point
                                % point k was shifted back by L/4
                                % interpolate between k and k-1
                                
                                if k ~= 1
                                    
                                    w_old = w_q_new(k-1,:) + (w_q_new(k,:) - w_q_new(k-1,:))./...
                                        (L_node(k)/2 + L_node(k-1)/2) * (L_node(k-1)/2 + 0.25 * L_node(k));
                                    
                                    g_old = g_q_new(k-1,:) + (g_q_new(k,:) - g_q_new(k-1,:))/...
                                        (L_node(k)/2 + L_node(k-1)/2) * (L_node(k-1)/2 + 0.25 * L_node(k));
                                    
                                else
                                    w_old = w_q_new(1,:);
                                    g_old = g_q_new(1,:);
                                end
                                
                                ind_div = k;
                                
                            else
                                % divide k+1 in half
                                
                                % values at the new point
                                w_star = w_q_new(k,:) + (w_q_new(k+1,:) - w_q_new(k,:))./...
                                    (L_node(k+1)/2 + L_node(k)/2) * (L_node(k)/2 + 0.25*L_node(k+1));
                                
                                g_star = g_q_new(k,:) + (g_q_new(k+1,:) - g_q_new(k,:))/...
                                    (L_node(k+1)/2 + L_node(k)/2) * (L_node(k)/2 + 0.25*L_node(k+1));
                                
                                % values at shifted k+1 point
                                % point k+1 was shifted up by L/4
                                % interpolate between k+1 and k+2
                                w_old = w_q_new(k+1,:) + (w_q_new(k+2,:) - w_q_new(k+1,:))./...
                                    (L_node(k+2)/2 + L_node(k+1)/2) * L_node(k+1)*0.25;
                                
                                g_old = g_q_new(k+1,:) + (g_q_new(k+2,:) - g_q_new(k+1,:))/...
                                    (L_node(k+2)/2 + L_node(k+1)/2) * L_node(k+1)*0.25;
                                
                                ind_div = k+1;
                                
                            end
                            
                            
                            % add to L_node, V_node
                            L_node(ind_div) = L_node(ind_div)/2;
                            L_node = [L_node(1:ind_div); L_node(ind_div); L_node(ind_div+1:end)];
                            
                            V_node(ind_div) = V_node(ind_div)/2;
                            V_node = [V_node(1:ind_div); V_node(ind_div); V_node(ind_div+1:end)];
                            
                            N_nodes = length(L_node);
                            
                            constants.L_node = L_node;
                            constants.V_node = V_node;
                            constants.N_nodes = N_nodes;
                            % now need to update y...
                            % and constants...
                            % and N_nodes (outside of constants)
                            
                            % find index of the kth node's w and g
                            
                            % have to put this in the right place
                            % for i = 1, it should be 1, 2, ... N_ab
                            % for i = 2, N_ab+1, N+2, ... N_ab+N_ab
                            % for i = i, it's (i-1)*N_ab, ... i*N_ab;%
                            ind_w_star = (k-1)*N_ab + [1:N_ab];
                            
                            % now have to add the temps and m, and U and
                            % such
                            ind_w_star = ind_w_star + 4+2*N_rw;
                            ind_g_star = (N_nodes-1)*N_ab + (k-1)*N_ab + [1:N_ab] + 4+2*N_rw;
                            
                            % index for the old point that's been shifted
                            % have to use (N_nodes-1) because N_nodes has
                            % already been incremenented
                            ind_w_old = (ind_div-1)*N_ab + [1:N_ab] + 4+2*N_rw;
                            ind_g_old = (N_nodes-1)*N_ab + (ind_div-1)*N_ab + [1:N_ab] + 4+2*N_rw;
                            
                            % put these into y_current as opposed to y_new
                            % because I'll have to redo the current step
                            
                            % insert the shifted existing point
                            y_current(ind_w_old) = w_old(:);
                            y_current(ind_g_old) = g_old(:);
                            
                            % insert the new point after the kth point
                            y_current = [y_current(1:ind_w_star(end)); ...
                                w_star(:); y_current(ind_w_star(end)+1:ind_g_star(end));...
                                g_star(:); y_current(ind_g_star(end)+1:end)];
                            
                            
                        end
                        fprintf('refined %0.d points: ', length(points_to_refine));
                        fprintf('%0.d, ', points_to_refine);
                        fprintf('N_nodes = %0.d\n', N_nodes);
                        
                        x_node(1) = L_node(1)/2;
                        for l = 2:N_nodes
                            x_node(l) = x_node(l-1) + (L_node(l) + L_node(l-1))/2;
                        end
                        
                        figure(48)
                        subplot(3,1,1)
                        plot(log10(L_node))
                        subplot(3,1,2)
                        plot(log10(w_q_new))
                        subplot(3,1,3)
                        plot(log10(g_q_new))
                    end
                    
					%  remove mesh points if error is too small in order to coarsen mesh
                    
                    %                 if length(points_to_coarsen) > 1
                    %                     % remove the zero
                    %                     points_to_coarsen = points_to_coarsen(2:end);
                    %                     spatial_error_OK = 0;
                    %                     error_OK = 0;
                    %                     n_coarsened = n;
                    %
                    %                     % use old points to refine space, not new points
                    %                     variables = unpack_y(y_current, constants);
                    %
                    %                     g_q_new = variables.g_q;
                    %                     w_q_new = variables.w_q;
                    %
                    %                     y_old = y_current;
                    %                     L_node_old = L_node;
                    %                     V_node_old = V_node;
                    %                     N_nodes_old = N_nodes;
                    %
                    %                     for j = 1:length(points_to_coarsen)
                    %
                    %                         % the point to coarsen (merge k and k+1)
                    %                         % need to correct for points already removed -> j-1
                    %                         k = points_to_coarsen(j) - (j - 1);
                    %
                    %                         % calculate values at new point
                    %                         w_bar(1,:) = w_q_new(k,:) + (w_q_new(k+1,:) - w_q_new(k,:))./...
                    %                             (L_node(k+1) + L_node(k)) * L_node(k+1);
                    %
                    %                         g_bar(1,:) = g_q_new(k,:) + (g_q_new(k+1,:) - g_q_new(k,:))/...
                    %                             (L_node(k+1) + L_node(k)) * L_node(k+1);
                    %
                    %
                    %                         % modify L_node, V_node
                    %                         L_node(k) = L_node(k) + L_node(k+1);
                    %                         L_node(k+1) = [];
                    %
                    %                         V_node(k) = V_node(k) + V_node(k+1);
                    %                         V_node(k+1) = [];
                    %
                    %                         N_nodes = length(L_node);
                    %
                    %                         constants.L_node = L_node;
                    %                         constants.V_node = V_node;
                    %                         constants.N_nodes = N_nodes;
                    %                         % now need to update y...
                    %                         % and constants...
                    %                         % and N_nodes (outside of constants)
                    %
                    %                         % find index of the kth node's w and g
                    %
                    %                         % have to put this in the right place
                    %                         % for i = 1, it should be 1, 2, ... N_ab
                    %                         % for i = 2, N_ab+1, N+2, ... N_ab+N_ab
                    %                         % for i = i, it's (i-1)*N_ab, ... i*N_ab;%
                    %                         ind_w_star = (k-1)*N_ab + [1:N_ab];
                    %
                    %                         % now have to add the temps and m, and U and
                    %                         % such, but remember N_nodes has been changed!
                    %                         ind_w_star = ind_w_star + 4+2*N_rw;
                    %                         ind_g_star = (N_nodes+1)*N_ab + (k-1)*N_ab + [1:N_ab] + 4+2*N_rw;
                    %
                    %                         % index for the old point that's been shifted
                    %                         % have to use (N_nodes-1) because N_nodes has
                    %                         % already been incremenented
                    %                         ind_w_kp1 = (k)*N_ab + [1:N_ab] + 4+2*N_rw;
                    %                         ind_g_kp1 = (N_nodes+1)*N_ab + (k)*N_ab + [1:N_ab] + 4+2*N_rw;
                    %
                    %                         % put these into y_current as opposed to y_new
                    %                         % because I'll have to redo the current step
                    %
                    %                         % modify the kth point
                    %                         y_current(ind_w_star) = w_bar(:);
                    %                         y_current(ind_g_star) = g_bar(:);
                    %
                    %                         % remove k+1 point
                    %                         y_current(ind_w_kp1) = [];
                    %                         y_current(ind_g_kp1) = [];
                    %
                    %                     end
                    %                     fprintf('removed %0.d points: ', length(points_to_coarsen));
                    %                 fprintf('%0.d, ', points_to_coarsen);
                    %                 fprintf('N_nodes = %0.d\n', N_nodes);
                    %
                    %                 end
                    
                end
                
            end
        else
            % not using adaptive scheme, don't need to check error
            error_OK = 1;
        end
        
    end
    
    
    
	% compute derivative at new value of y
	% also extract all the variables I want to store that are normally only found within the differential equations
    
	constants.step = 1;
    
    try
        
        [f_np1, debug_data] = diffeqns(y_new, constants, guesses, PDT);
        
    catch
        
        f_np1 = f;
        f_np1(1) = nan;
        disp('error in calling diffeqns to get the debug data, not integrating')
    end
    
    diff_eqns_error_flag = 0;
    
    m_tg(n+1) = debug_data.m_tg;
    U_tg(n+1) = debug_data.U_tg;
    m_l(n+1) = debug_data.m_l;
    T_l(n+1) = debug_data.T_l;
    g_q{n+1} = debug_data.g_q;
    r_q{n+1} = debug_data.r_q;
    mom{n+1} = debug_data.mom;
    w_q{n+1} = debug_data.w_q;
    V_bubi{n+1} = debug_data.V_bubi;
    T_s(n+1) = debug_data.T_s;
    P(n+1) = debug_data.P;
    deltaT_sup(n+1) = debug_data.deltaT_sup;
    rho_l(n+1) = debug_data.rho_l;
    rho_tg(n+1) = debug_data.rho_tg;
    rho_tg_sat(n+1) = debug_data.rho_tg_sat;
    x_tg(n+1) = debug_data.x_tg;
    V_l(n+1) = debug_data.V_l;
    V_tg(n+1) = debug_data.V_tg;
    V_bub(n+1) = debug_data.V_bub;
    m_bub(n+1) = debug_data.m_bub;
    U_bub(n+1) = debug_data.U_bub;
    V_l_star(n+1) = debug_data.V_l_star;
    n_bubi(n+1) = debug_data.n_bubi;
    fill_level(n+1) = debug_data.fill_level;
    Qdot_lw(n+1) = debug_data.Qdot_lw;
    gas_holdup(n+1) = debug_data.gas_holdup;
    N_full(n+1) = debug_data.N_full;
    U_liq(n+1) = debug_data.U_liq;
    mdot_out_liq(n+1) = debug_data.mdot_out_liq;
    mdot_out_vap(n+1) = debug_data.mdot_out_vap;
    h_l(n+1) = debug_data.h_l;
    h_tg_sat(n+1) = debug_data.h_tg_sat;
    gas_holdup_injector(n+1) = debug_data.gas_holdup_injector;
    Vdot_bub(n+1) = debug_data.Vdot_bub;
    T_lw(:,n+1) = debug_data.T_lw;
    T_gw(:,n+1) = debug_data.T_gw;
    u_rise{n+1} = debug_data.u_rise;
    N_bubi{n+1} = debug_data.N_bubi;
    mdot_tg(n+1) = debug_data.mdot_tg;
    mdot_bub_l{n+1} = debug_data.mdot_bub_l;
    mdot_bub_tg{n+1} = debug_data.mdot_bub_tg;
    mdot_l(n+1) = debug_data.mdot_l;
    Cp_l(n+1) = debug_data.Cp_l;
    node_level{n+1} = debug_data.node_level;
    diff_eqns_error_flag = debug_data.diff_eqns_error_flag;
    
    
    
    %     fprintf('r_q = ')
    %     fprintf('%4.6g,\t', r_q(1:6,1))
    %     fprintf('\n')
    %
    %     fprintf('w_q = ')
    %     fprintf('%4.6g,\t ', w_q(1:6,1))
    %     fprintf('\n')
    
    %     CFL = h*max(u_rise(:))/L_node;
    %
    %     if CFL > 1
    %         disp('CFL > 1')
    %     end
    %
    if max(r_q{n+1}(:)) > 0.99
        disp('r_q is getting too big')
    end
    
    
    if strcmp(ADQMOM, 'on')
        
        
        % possible values of p
        p_max = (N_mom - 1);
        p_vec = [1:p_max]/3;
        
        for i = 1:length(p_vec)
            
            % equivalent abscissas
            l = r_q.^(1/p_vec(i));
            l_norm = l/max(l);
            
            l_std(i) = std(l_norm);
            
        end
        
        % take p that has max std
        [~, i_p] = max(l_std);
        
        p = p_vec(i_p);
        
        disp(num2str(p))
        
        constants.ADQMOM_p = p;
        V_moment_index = 1 + 3*p;
        constants.V_moment_index = V_moment_index;
        
    elseif strcmp(ADQMOM,'off')
        p = 1;
    end
    
    
    if plot_periodically
        
        clock_now = clock;
        
        time_since_plot = etime(clock_now, clock_plot);
        
        if time_since_plot >= t_plot;
            
            figure(87)
            clf
            
            % number density, gas holdup, sauter mean diameter
            x_node = [1:length(V_bubi{end})];
            plot(x_node, N_bubi{end}, 'k-s', x_node, V_bubi{end}, 'k-o', x_node, (6*V_bubi{end}./N_bubi{end}/pi).^(1/3),'k-*')
            set(gca,'yscale','log')
            legend('number density','volume density','SMD')
            
            
            clock_plot = clock;
        end
    end
    
    
    
    guesses.P = P(n+1);
    guesses.rho_tg = rho_tg(n+1);
    guesses.rho_l = rho_l(n+1);
    guesses.node_level = node_level{n+1};
    
    t(n+1) = t(n) + h;
    
    n = n + 1;
    
    y_current = y_new;
    
    if deltaT_sup(n) > 7.5
        running = 0;
        stop_reason = 'superheat got bigger than 7.5';
    end
    
    if (sum(abs(imag(y_new))) > 0) || (sum(isnan(y_new)) > 0)
        running = 0;
        disp('imaginary or nans')
        stop_reason = 'imaginary or nans in y(:,n)';
    end
    
    if (t(n) > t_end) || ( fill_level(n) < 1e-6)
        %         disp('reached end t or ran out of liquid')
        running = 0;
        stop_reason = 'reached t_end or ran out of liquid';
    end
    
    
    if P(end) < 2e5
        running = 0;
        %         disp('pressure got low')
        stop_reason = 'pressure got below ambient';
    end
    
    gas_holdup_next = gas_holdup(end) + (gas_holdup(end) - gas_holdup(end-1));
    
    if gas_holdup_next > 0.90
        running = 0;
        %         disp('alpha -> 0.8')
        stop_reason = 'alpha got big';
    end
    
    clock_now = clock;
        
    time_since_save = etime(clock_now, clock_save);
    
    if time_since_save >= t_save && (save_periodically && save_stuff)
        
        current_dir = pwd;
        cd(results_save_dir)
        
        clock_save = clock;
        if first_save == 1
            
            file_exists = 1;
            file_num = 1;
            
            while file_exists
                
                save_filename_intermediate_numbered = [save_filename num2str(file_num) '_intermediate_data.mat'];
                save_filename_numbered = [save_filename num2str(file_num) '.mat'];
                
                if ~exist(save_filename_numbered,'file') && ~exist(save_filename_intermediate_numbered,'file')
                    disp('saving intermediate data')
                    save(save_filename_intermediate_numbered,'-v7.3')
                    file_exists = 0;
                    first_save = 0;
                    
                end
                
                file_num = file_num+1;
                
            end
        
        else
            
            disp('saving')
            save([save_filename_intermediate_numbered],'-v7.3')
        end
        
        cd(current_dir)
            
        
    end
    
    if time_out
        clock_now = clock;
        time_since_start = etime(clock_now, clock_start);
        if time_since_start > max_comp_time
            running = 0;
            stop_reason = 'timed out';
        end
    end
    
    if constants.error_detected
        runnung = 0;
        stop_reason = 'error detected';
    end
      
end
if exist('stop_reason', 'var')
    disp(stop_reason)
else
    disp('no stop reason assigned')
end

if h_min_error_count > 0
    disp(['exceeded min step size ' num2str(h_min_error_count) ' times'])
end

%% plotting and output

if save_stuff == 1
    % store the current directory so we can get back to it later
    current_dir = pwd;
    
    cd(results_save_dir)
    
    if save_parameters_only
        
        if ~exist('n_min','var')
            n_min = 1;
        end
        
        if ~exist('n_peak','var')
            n_peak = 1;
        end
        t_min = t(n_min);
        t_peak = t(n_peak);
        t_LRO = t(end);
        P_min = P(n_min);
        P_peak = P(n_peak);
        P_LRO = P(end);
        alpha_f = V_bubi(end);
        fill_level_f = fill_level(end);
        P_lin = (P_LRO - P_peak)/(t_LRO - t_peak)*(t - t_peak) + P_peak;
        P_dev = mean( abs(P(n_peak:end) - P_lin(n_peak:end)) );
        
        save(save_filename,'t_min','t_peak','t_LRO'...
            ,'P_min','P_peak','P_LRO'...
            ,'alpha_f','fill_level_f','P_dev','-v7.3')
        
    else
        
        if exist('save_filename_intermediate_numbered','var')
            save(save_filename_numbered,'-v7.3')
        else
        
            file_exists = 1;
            file_num = 1;
            while file_exists
                save_filename_numbered = [save_filename num2str(file_num) '.mat'];
                if ~exist(save_filename_numbered,'file')
                    save(save_filename_numbered,'-v7.3')
                    file_exists = 0;
                end
                file_num = file_num+1;
            end
        end
        
    end
    
    cd(current_dir);
    
end



if plot_stuff == 1
    
    m_out = cumtrapz(t, mdot_out_liq + mdot_out_vap);
    mh_out = cumtrapz(t, (mdot_out_liq.*h_l + mdot_out_vap.*h_tg_sat));
    
    [P_exp, T_lw_out_exp, LL_exp] = load_experimental_data(t, specified_case);
    
    %     T_l = y(4+N_rw,:);
    %     m_l = y(3+N_rw,:);
    
    %     m_tg = y(1,:);
    
    T_lw_in = T_lw(1,:);
    T_gw_in = T_gw(1,:);
    T_lw_out = T_lw(end,:);
    
    %         m_tg = y(1);
    % U_tg = y(2);
    % T_gw = y(3 : 2+N_rw);
    % m_l = y(3+N_rw);
    % T_l = y(4+N_rw);
    % T_lw = y(5+N_rw : 4+2*N_rw);
    
    
    %
    figure(1)
    hold on
    plot(t,P/1e6,'k-')
    plot(t,P_exp/1e6,'k--')
    legend('model','experiment')
    xlabel('Time [s]')
    ylabel('Pressure [MPa]')
    %         hold on
    %         if exist('t_peak','var')
    %             plot(t_min, P(n_min)/1e6, 'ko')
    %             plot(t_peak, P(n_peak)/1e6, 'ks')
    %         end
    
    title('pressure')
    
    figure(2)
    hold on
    plot(t,T_l,'k-',t,T_s,'r:')
    legend('Liquid','T_{sat}(P) = T_{tg}')
    ylabel('Temperature')
    xlabel('Time [s]')
    title('temperatures')
    
    figure(3)
    hold on
    plot(t,T_lw_in,'k-',t,T_gw_in,'b--',t,T_lw_out,'r-.',t,T_lw_out_exp+273.15,'k:')
    title('wall temp')
    xlabel('Time [s]')
    legend('liquid','vapor','outside liquid','experimental outside liquid')
    title('wall temp')
    
    figure(4)
    hold on
    plot(t,m_l,'k-',t,m_tg,'b--',t,m_bub,'k:',t, m_out ,'r--', t, m_l + m_tg + m_bub + m_out, 'g--')
    title('Mass')
    xlabel('Time [s]')
    legend('Liquid','Vapor','Bubbles','Out through injector','Sum')
    title('masses')
    
    figure(5)
    hold on
    plot(t, x_tg)
    xlabel('Time [s]')
    ylabel('vapor mass fraction []')
    title('ullage vapor mass fraction')
    
    
    figure(6)
    hold on
    plot(t, fill_level, 'k', t, LL_exp, 'b')
    xlabel('Time [s]')
    ylabel('fill level [%]')
    title('fill level')
    legend('Model','Experiment')
    
    
    %         figure(7)
    %         hold on
    %         plot(t, A_bub,'k')
    %         xlabel('Time [s]')
    %         ylabel('A/V [1/m]')
    %         title('interfacial area per volume')
    
    figure(8)
    hold on
    plot(t, V_bub./V_l_star, 'k', t, gas_holdup_injector)
    xlabel('Time [s]')
    ylabel('gas holdup')
    title('gas holdup')
    legend('mean','bottom of tank')
    
    %         figure(9)
    %         hold on
    %         plot(t, 6*V_bub./A_bub, 'k')
    %         xlabel('Time [s]')
    %         ylabel('sauter mean diameter [m]')
    %         set(gca,'yscale','log')
    %         title('sauter mean diameter')
    
    %         figure(10)
    %         hold on
    %         plot(t, r_q)
    %         xlabel('Time [s]')
    %         ylabel('abscissas [m]')
    %         set(gca,'yscale','log')
    %         title('abscissas')
    
    %         figure(11)
    %         hold on
    %         plot(t, w_q)
    %         xlabel('Time [s]')
    %         ylabel('weights [?]')
    %         set(gca,'yscale','log')
    %         title('weights')
    
    figure(12)
    hold on
    plot(t, V_l, 'k', t, V_tg, 'k:', t, V_bub, 'k--', t, V_l+V_tg+V_bub,'k-.')
    xlabel('Time [s]')
    ylabel('volume [m^3]')
    title('volumes')
    legend('liquid (pure)','ullage','bubbles','sum')
    
    figure(13)
    hold on
    plot(t, n_bubi, 'k')
    xlabel('Time [s]')
    ylabel('number/m^3')
    title('bubble number density')
    
    figure(14)
    hold on
    plot(t, mh_out, 'k', t, U_liq, 'k:', t, U_tg, 'k--', t, U_bub, 'k-.', t, (mh_out + U_liq + U_tg + U_bub), 'b')
    xlabel('Time [s]')
    ylabel('energy [J]')
    legend('out the injector','liquid','vapor','bubbles','sum')
    
    figure(15)
    hold on
    plot(t, [0; diff(t(:))], 'k')
    xlabel('Time [s]')
    ylabel('dt [s]')
    title('time step')
    set(gca,'yscale','log')
    
    figure(16)
    hold on
    plot(t, ind_max_rel_err(1:length(t)))
    xlabel('Time [s]')
    ylabel('index []')
    title('index of max relative error')
    
    
    
end

beep

toc

%% differential equations
function varargout = diffeqns(y, constants, guesses, PDT)

% retrieve constants
D_tank = constants.D_tank;
t_w = constants.t_w;
rho_w = constants.rho_w;
cv_w = constants.cv_w;
Cd = constants.Cd;
A_inj = constants.A_inj;
T_air = constants.T_air;
V_tank = constants.V_tank;
k_w = constants.k_w;
K_b = constants.K_b;
N_A = constants.N_A;
P_cr = constants.P_cr;
T_cr = constants.T_cr;
h_planck = constants.h_planck;
g = constants.g;
C_hamaker = constants.C_hamaker;

n_nuc_freq = constants.n_nuc_freq;

phi = constants.phi;
C_rdot = constants.C_rdot;
C_nuc_rate = constants.C_nuc_rate;
C_r_nuc = constants.C_r_nuc;
C_coalescence = constants.C_coalescence;
C_dTs = constants.C_dTs;
C_qdot_lw = constants.C_qdot_lw;

L_node = constants.L_node;
% V_node = 0.25*pi*D_tank^2*L_node;
V_node = constants.V_node;
N_nodes = constants.N_nodes;
fluid = constants.fluid;
N_ab = constants.N_ab; % number of abscissas
N_mom = 2*N_ab;
hesson_fit = constants.hesson_fit;
% bubble_rise_velocity_fit = constants.bubble_rise_velocity_fit;

p = constants.ADQMOM_p;	% parameter for ADQMOM (not really supported right now)
Ru = 8314.4;	% universal gas constant

% extract the velocity of the liquid level
u_LL = guesses.dLL_dt;

% check for nan's and set to 0 if it is nan
if isnan(u_LL)
    u_LL = 0;
end

% retrieve derivatives calculated with backwards differencing
% not currently in use
% Pdot = derivatives(1);
% rhodot_l = derivatives(2);
% rhodot_tg = derivatives(3);
% Vdot_l = derivatives(4);
% Vdot_tg = derivatives(5);

% check for negative values
% if we find any, just multiply them by -1, store the indices, and print something to screen
if min(y) < 0
    ind_negative = find(y < 0);
    y(ind_negative) = -y(ind_negative);
    fprintf('negative part of y. # = %0.d, ',length(ind_negative))
    if length(ind_negative > 1)
        fprintf('indices: ')
        fprintf('%0.d, ', ind_negative)
        fprintf('\n')
    else
        fprintf('index: %0.d\n',ind_negative)
    end
    
else
    ind_negative = [];
end

% unpack all the variables from the "y" vector
variables = unpack_y(y, constants);

m_tg = variables.m_tg;	% ullage mass
U_tg = variables.U_tg;	% ullage internal energy (not specific)
T_gw = variables.T_gw;  % temperatures in the wall for ullage
m_l = variables.m_l;	% mass of liquid
T_l = variables.T_l;	% liquid temperature
T_lw = variables.T_lw;	% temperatures in the wall for liquid

g_q = variables.g_q;	% weighted abscissas
w_q = variables.w_q;	% weights

% error checking for imaginary or nans
if isnan(sum(y)) || ~isreal(sum(y))
    disp('problem: nans or imaginary y')
end

% get the abscissas from the weighted abscissas
r_q = g_q./w_q;

% error check for negative abscissas
if min(r_q) < 0
    disp('negative abscissa')
end

% get the moments from gauss quadrature
for i = 1:2*N_ab
    mom(:, i) = sum( r_q.^((i-1)/p) .* w_q, 2 );
end

mom(mom<0) = 0; % if moments are negative anywhere, set to 0

N_bubi = mom(:,1); % the number density of bubbles (0th moment)

% error check for imaginary moments
if ~isreal(mom)
    disp('imaginary moments')
end

% index of the moment for r^3 (what is used to calculate volume)
V_moment_index = constants.V_moment_index;

% bubble volume per unit volume of liquid/bubble mixture (hence the i)
% (can also view this as the vapor volume fraction aka gas holdup aka void fraction)
V_bubi = 4/3*pi*mom(:,V_moment_index);

% error checking for imaginary values of the void fraction
if sum(imag(V_bubi)) > 0
    disp('V_bubi went imaginary')
end

% see if the void fraction is greater than 1 anywhere
% if so, set an error flag and set the 
% void fraction to 0.99
if max(V_bubi) > 1
    ind_V_bubi_gt1 = find(V_bubi > 1);
    V_bubi(V_bubi>1) = 0.99;
    V_bubi_error_flag = 1;
else
    V_bubi_error_flag = 0;
end

% get system pressure
% (assumes pressure is same throughout tank, with no gravity head)
P = get_P_from_mU_mT(m_tg, U_tg, m_l, T_l, V_bubi, PDT, constants, guesses);

% check to see if there was an error in the function that calculated P
% (pi is an error flag)
if (P == pi) || isnan(P)
    disp('P error')
    constants.error_detected = 1;
    P = guesses.P;
end


% error checking for supercritical pressure
if P > constants.P_cr
    disp('P greater than critical')
    constants.error_detected = 1;
    P = guesses.P;
end

% get density of liquid and ullage based on temperature and pressure
% get it either from interpolation of a table or by calling refprop

if strcmp(constants.property_source,'PDT')
    % interpolating a table
	
	% liquid density
    rho_l = qinterp2(PDT.T, PDT.P, PDT.D_liq, T_l, P/1e3);
    
    % if rho_l is NaN, it means we went outside the bounds of PDT, so
    % instead extrapolate it using interp2 (slower than qinterp2)
    if isnan(rho_l)
        rho_l = interp2(PDT.T, PDT.P, PDT.D_liq, T_l, P/1e3, 'spline');
    end
    
	% get other properties from a fit based on P
    [rho_tg_l, rho_tg_v, u_tg_l, u_tg_v] = fits_for_getting_P(P, fluid);
    
    
elseif strcmp(constants.property_source,'refprop')
	% calling refprop
    
	% liquid density (have to solve nonlinear equation for it)
    rho_l = get_D_from_TP(T_l, P, guesses, constants, fluid);
    
	% get other properties from refprop
    [rho_tg_l, rho_tg_v, u_tg_v] = refpropm('+-U','P',P/1e3,'Q',1,fluid);
    u_tg_l = refpropm('U','P',P/1e3,'Q',0,fluid);
    
end


% get saturation properties for ullage

% internal energy (specific)
u_tg = U_tg/m_tg;

% vapor mass fraction
x_tg = (u_tg - u_tg_l)/(u_tg_v - u_tg_l);

% void fraction (vapor volume fraction)
alpha = 1/( 1 + rho_tg_v/rho_tg_l * (1 - x_tg)/x_tg );

% density (mixture)
rho_tg = alpha*rho_tg_v + (1 - alpha)*rho_tg_l;

% just used for debugging
% if constants.step == 1
%     fprintf('1 - x_tg = %0.3g\n',1 - x_tg)
% end

% error checking for x_tg > 1 (means that it's now superheated vapor)
% if x_tg > 1
%     disp('x_tg > 1')
%     x_tg = 1;
% end

% total liquid volume (not including bubbles)
V_l = m_l./rho_l;

% how full each node is. will be 1 for all nodes below liquid level, 0 for
% all nodes above, and [0,1] for the node in which the liquid level
% currently sits
node_level = get_node_levels(V_l, V_bubi, V_node);%, guesses.node_level);

% total bubble volume
V_bub = sum_over_nodes(V_bubi, node_level, V_node);

% total volume of liquid + bubbles
V_l_star = V_l + V_bub;

% for reference:
% V_l_star = volume of liquid and bubbles
% = V_l + V_bub
% V_bub = V_bubi*(V_l + V_bub)
% V_l_star = V_l./(1 - V_bubi);

% number of full nodes still in the liquid
N_full = sum(node_level == 1);

% check if an error was thrown earlier b/c void fraction > 1
% if it happened at a node with no liquid in it, set it to 0
% otherwise, leave as is (it was set to 0.99 earlier)
if V_bubi_error_flag
    if min(ind_V_bubi_gt1) < N_full + 1
        %V_bubi is > 1 at a node we care about
    else
        V_bubi_error_flag = 0;
    end
end

% volume of ullage
V_tg = m_tg/rho_tg;

% liquid properties

% metastable liquid properties (at current P and T_l)
[h_l, dh_drho_l, drho_dP_l, u_l, Cv_l, dP_dT_l, ...
    k_l, Cp_l, s_l, MW, mu_l] = ...
    refpropm('H!RUO#LCSMV','T',T_l,'D&',rho_l,fluid);


% get properties of liquid saturated at current liquid temperature
[P_sat, s_liq_sat, h_liq_sat] = refpropm('PSH', 'T', T_l, 'Q', 0, fluid);

dP_drho_l = 1e3./drho_dP_l; % partial derivative of P wrt rho (don't know what's being held constant -> look at refprop.m)
dP_dT_l = dP_dT_l*1e3; % partial derivative of P wrt T (not sure what's being held constant -> look at refprop.m)
alpha_l = k_l./(rho_l .* Cp_l); % thermal diffusivity
P_sat = 1e3*P_sat; % saturation pressure (kPa->Pa)

nu_l = mu_l/rho_l; % liquid viscosity (dynamic I think, although I always get dynamic/kinematic confused)

Pr_l = nu_l/alpha_l; % Prandtl number of the liquid


% need to calculate a bunch of partial derivatives for the Vdot function.
% to do so, need to get a bunch of stuff out of refprop. To sove time I also
% get other properties that will be needed elsewhere (refprop = slow)

% refprop is first called to get saturated liquid properties, then for saturated vapor
% (saturated at the current pressure)

% first call
% properties needed for Vdot calculation (liquid)
[u_tg_l_sat, rho_tg_l, dP_dT_tg_sat, drho_dP_T, drho_dT_P, dh_dT_P, dh_dP_T...
    ,sigma, h_l_sat, T_tg] = refpropm('UDERW(*IHT', 'P', P/1e3, 'Q', 0, fluid);
	
% partial derivative of P wrt T along saturation (convert kPa to Pa)
dP_dT_tg_sat = dP_dT_tg_sat * 1e3;

% partial derivative of rho wrt P at constant T (convert kPa to Pa)
drho_dP_T = drho_dP_T * 1e-3;

% partial derivative of h wrt P at constant T (convert kPa to Pa)
dh_dP_T = dh_dP_T * 1e-3;

% partial derivative of u wrt T at constant P
du_dT_P = dh_dT_P + P/rho_tg_l^2 * drho_dT_P;

% partial derivative of u wrt P at constant T
du_dP_T = dh_dP_T + 1/rho_tg_l + P/rho_tg_l^2 * drho_dP_T;

% partial derivative of u wrt T along saturation (liquid)
du_dT_sat_tg_l = du_dT_P + du_dP_T * dP_dT_tg_sat;

% partial derivative of rho wrt T along saturation (liquid)
drho_dT_l_sat = drho_dT_P + drho_dP_T * dP_dT_tg_sat;

% partial derivative of rho wrt P along saturation (liquid)
% drho_dP_sat = drho_dP_T + drho_dT_P / dP_dT_sat; % don't need this anymore with new Vdot fn


% second call
% properties needed for Vdot calculation (vapor)
[u_tg_v_sat, rho_tg_v, drho_dP_T, drho_dT_P, dh_dT_P, dh_dP_T,...
    T_s, h_tg_sat, s_tg_sat, Cp_tg, mu_tg] = ...
    refpropm('UDRW(*THSCV', 'P', P/1e3, 'Q', 1, fluid);
	
% partial rho / partial P at constant T (convert units, refprop does kPa)
drho_dP_T = drho_dP_T * 1e-3;

% partial h / partial P at constant T (convert units, refprop does kPa)
dh_dP_T = dh_dP_T * 1e-3;

% partial u / partial T at constant P
du_dT_P = dh_dT_P + P/rho_tg_v^2 * drho_dT_P;

% partial u / partial P at constant T
du_dP_T = dh_dP_T + 1/rho_tg_v + P/rho_tg_v^2 * drho_dP_T;

% partial u / partial T along saturation line (vapor)
du_dT_sat_tg_v = du_dT_P + du_dP_T * dP_dT_tg_sat;

% partial rho / partial T along saturation line (vapor)
drho_dT_v_sat = drho_dT_P + drho_dP_T * dP_dT_tg_sat;

% partial rho / partial P along saturation line (vapor)
% drho_dP_sat_tg = drho_dP_T + drho_dT_P / dP_dT_sat; % don't need this anymore with new Vdot fn

% partial rho / partial x at constant P
drho_dx_P_tg = -rho_tg^2 *(1/rho_tg_v - 1/rho_tg_l);

% partial rho / partial P at constant x
drho_dP_x_tg = (1/dP_dT_tg_sat) * rho_tg^2 * ( x_tg/rho_tg_v^2 * drho_dT_v_sat + ...
    (1-x_tg)/rho_tg_l^2 * drho_dT_l_sat );

rho_tg_sat = rho_tg_v; % density of saturated vapor (single phase)

% heat of vaporization (at saturation based on P)
h_lv = h_tg_sat - h_l_sat;

% liquid superheat
deltaT_sup = T_l - T_s;

% vapor mass fration going out the injector
x_out = constants.C_x_inj * V_bubi(1)/(V_bubi(1) + rho_l/rho_tg_sat*(1 - V_bubi(1)));

% density, entropy, and enthalpy of the flow going out the injector
rho_liq_mix = rho_l*(1-V_bubi(1)) + rho_tg_sat*V_bubi(1);
s_liq_mix = s_l * (1-x_out) + s_tg_sat * x_out;
h_liq_mix = h_l * (1-x_out) + h_tg_sat * x_out;

% mass flow rate out through the injector
mdot_out_mix = A_inj*Cd*injector_flow(x_out, P, hesson_fit, rho_liq_mix, h_liq_mix, fluid, constants);

% the liquid and vapor flow rates
mdot_out_liq = (1 - x_out)*mdot_out_mix;
mdot_out_vap = x_out*mdot_out_mix;

% bulk flow velocity out the bottom
% I was using a constant that I set as an input parameter, but I commented that out when I started picking and choosing various terms to include u_bulk, not just as an all-or-nothing thing
% if constants.include_u_bulk
u_bulk = mdot_out_mix / rho_l / (0.25 * pi * D_tank^2);
% else
%     u_bulk = 0;
% end

% Morton number
Mo = g*mu_l^4*(rho_l - rho_tg_sat)/(rho_l^2*sigma^3);

% parameters for fan-tsuchiya rise velocity expression
% (not in use currently)
% n_FT = 1.6;
% c_FT = 1.2;
% Kbo_FT = 10.2;
% Kb_FT = max([ Kbo_FT*Mo^-0.038, 12]);

% initialize u_rise
u_rise = zeros(size(r_q));

% difference between liquid and vapor density (used frequently)
delta_rho = rho_l - rho_tg_sat;

% loop through the nodes and calculate rise velocity using ishii & zuber, 1979
for i = 1:N_full+1
    
	% check if the void fraction is high enough to be churn-turbulent
    if V_bubi(i) > constants.CD_churn_turb_boundary
		% if it is, use this simple expression
        Eprime = (1 - V_bubi(i)).^2;
        Cd = 8*Eprime/3;
        u_rise(i,:) = sqrt(8/3 * delta_rho * g * r_q(i,:)./(rho_l * Cd));
    else
        % not churn turbulent based on void fraction (could still be based on other criteria)
		
		% need to check my thesis/the ishii and zuber paper to remember the regime names
        
        alpha_IZ = V_bubi(i)*ones(1,N_ab); % turn void fraction into vector spanning each abscissa/weight
        
		% mixture viscosity
        mu_mix = mu_l*(1 - alpha_IZ).^(-2.5*(mu_tg + 0.4*mu_l)/(mu_tg + mu_l));
        
        % a guess at u, found from assuming C_D = 24/Re
        % u_guess = 1/9 * r_q(i,:).^2 * (rho_l - rho_v) * g/mu_mix;
        
		% viscosity number used by ishii and zuber
        N_mu = Mo^(1/4);
        
		% non-dimensional radius
        r_d_star = r_q(i,:) * (rho_l * g * delta_rho/mu_l^2)^(1/3);
		
		% empirical function of the non-dim. radius
        psi = 0.55*( (1 + 0.08 * r_d_star.^3).^(4/7) - 1).^(3/4);
		
		% rise velocity
        u_rise(i,:) = 10.8*mu_l./(rho_l * r_q(i,:)) .* mu_l./mu_mix .*(1-alpha_IZ).^2 ...
            .* psi.^(4/3) .* (1 + psi) ./ ...
            (1 + psi.*( mu_l./mu_mix .*(1 - alpha_IZ).^0.5).^(6/7) );
        
        
		% alternative:
        % fit from fan and tsuchiya to get the single particle rise u
        % then use ishii and zuber to get the mixture result
        
        %         De = 2*r_q(i,:);
        %
        %         De_ND = De*sqrt(rho_l*g/sigma);
        %
        %
        %         u_rise_ND = ( (Mo^-0.25/Kb_FT * ( delta_rho/rho_l )^(5/4) * De_ND.^2 ).^-n_FT ...
        %             + ( 2*c_FT./De_ND + ( delta_rho/rho_l ) * De_ND/2 ).^(-n_FT/2) ).^(-1/n_FT);
        %
        %         u_rise_FT = u_rise_ND/(rho_l/(sigma*g))^(1/4);
        %
        %         % u_rise = constants.C_u_rise * u_rise;%.*(( (1 - V_bubi.^(5/3))./(1 - V_bubi).^2) * ones(size(u_rise(1,:))) ).^(-0.5);
        %
        %         alpha_IZ = V_bubi(i)*ones(1,N_ab);
        %
        %         mu_mix = mu_l*(1 - alpha_IZ).^(-2.5*(mu_tg + 0.4*mu_l)/(mu_tg + mu_l));
        %
        %         Re_inf = rho_l * u_rise_FT .* 2 .* r_q(i,:)./mu_l;
        %
        %         u_rise(i,:) = u_rise_FT.* (1 - alpha_IZ).*mu_l./mu_mix .* ...
        %             (1 + 0.1*Re_inf.^0.75)./(1 + 0.1*Re_inf.^0.75.*( sqrt(1 - alpha_IZ).*mu_l./mu_mix).^(6/7));
        %
        %
        %         N_mu = Mo^(1/4);
        
		
		% now loop through each weight/abscissa
        for j = 1:N_ab
            
			% non-dimensional radius
            r_d_star = r_q(i,j) * (rho_l * g * delta_rho/mu_l^2)^(1/3);
			
			% empirical function of the non-dim. radius
            psi = 0.55*( (1 + 0.08 * r_d_star^3)^(4/7) - 1)^(3/4);
            
            % compare psi to N_mu to see if we're in a different regime (if not, use already calculated value for u_rise)
            if N_mu > 0.11 * (1 + psi)/psi^(8/3)
                
				% function of the void fraction that's useful
                f_alpha = sqrt(1 - V_bubi(i)) * mu_l/mu_mix(j);
				
				% 
                E = ( (1 + 17.67*f_alpha.^(6/7) )./(18.67*f_alpha)).^2;
				
				% Eotvos number
                Eo = g*delta_rho*4*r_q(i,j).^2/sigma;
				
				% drag coefficient
                Cd = 2/3*E*sqrt(Eo);
                
				% Cd maxes out (spherical cap bubble I think)
                if Cd > 8/3  * (1 - V_bubi(i))^2
                    Cd = 8/3  * (1 - V_bubi(i))^2;
                end
                
				% rise velocity
                u_rise(i,j) = sqrt(8/3 * delta_rho * g * r_q(i,j)./(rho_l * Cd));
                
				% alternative method: use fzero to solve a different ishii and zuber expression
                %             else
                %                         u_guess = 1/9 * r_q(i,j).^2 * delta_rho * g/mu_mix(j);
                %
                %                         ishii_zuber_fn = @(u) 24/(rho_l*r_q(i,j)*u/mu_mix(j))*(1 + 0.1*(rho_l*r_q(i,j)*u/mu_mix(j))^0.75) - 8/3*r_q(i,j)*delta_rho*g/(u^2*rho_l);
                %                         u_rise(i,j) = fzero(@(u)ishii_zuber_fn(u), u_guess);
                %
            end
        end
        
    end
    
end


% set u_rise for nodes with nothing in them to 0
u_rise(N_full+2:end,:) = 0;

% add in a constant that multiplies u_rise (usually = 1)
u_rise = constants.C_u_rise*u_rise;

% subtract off the bulk velocity
% u_rise = u_rise - u_bulk;

% compute the CFL and check if it's too big
if constants.step == 1
    %     u_max = max(u_rise(:));
    CFL = max(constants.h*u_rise./(L_node*ones(1,N_ab)));
    if CFL > 0.5
        fprintf('max CFL>0.5, = %0.3g\n', CFL);
    end
end

% min and max Reynolds numbers (for debugging)
Re_max = max(max(r_q(1:N_full+1,:).*u_rise(1:N_full+1,:)))*rho_l/mu_l;
Re_min = min(min(r_q(1:N_full+1,:).*u_rise(1:N_full+1,:)))*rho_l/mu_l;


% Eotvos number and the min and max values in the tank (for debugging)
Eo = g*(rho_l - rho_tg_sat)*r_q.^2*4/sigma;

Eo_max = max(max(Eo(1:N_full+1,:)));
Eo_min = min(min(Eo(1:N_full+1,:)));


% this is handy for debugging
% if constants.step == 1 && constants.t > 0.001
%     fprintf('Re: %0.2g - %0.2g, Eo: %0.2g - %0.2g, Mo: %0.2g, Ja: %0.2g\n', Re_min, Re_max, Eo_min, Eo_max, Mo, Ja)
% end

% u_rise(r_q < 1e-6) = 0;

% some error checking
if sum(sum(isnan(u_rise(1:N_full+1,:)))) > 0
    disp('nans in u_rise')
end

% compute superficial gas velocity
% needed for computing the bubble diffusivity
for i = 1:N_full + 1
    u_superficial(i) = 4/3 * pi * sum(r_q(i,:).^(3) .* w_q(i,:) .* (u_rise(i,:)) );
end


% axial diffusivity (Hikita and Kikukawa, 1974)
% (zz is supposed to mean along the length of the tank, even though I use x elsewhere...)
D_zz = (0.15 + 0.69*abs(u_superficial).^0.77)*D_tank^1.25 * mu_l^-0.12;
% other expressions - page 185, #6


% initialize the spatial derivatives

duw_dx = zeros(size(w_q));
dug_dx = duw_dx;

% create some vectors that will be handy later on
u_vec = u_rise;
uw_vec = (u_rise).*w_q;
ug_vec = (u_rise).*g_q;

% u_vec = u_rise - u_bulk;
% uw_vec = (u_rise - u_bulk).*w_q;
% ug_vec = (u_rise - u_bulk).*g_q;

% fourier_number = mean(D_zz) * constants.h / mean(L_node)^2;


% compute spatial derivative terms
% 1st order derivative for convection
% 2nd order derivative for diffusion

for i = 1:N_full + 1
    
    % diffusion terms
    if N_full >= 1
		% if there's at least 1 full node (so I can use i + 1)
		
		% interior points
        if i > 1 && i < N_full + 1
            
			% diffusion constant (D) at i + 1/2
            D_ip2 = 0.5*( D_zz(i+1) + D_zz(i) );
			
			% now at i-1/2
            D_im2 = 0.5*( D_zz(i) + D_zz(i-1) );
			
			% and at i
            D_i = D_zz(i);
            
			% d/dx( D * (dw/dx) )
            dDdw_dx2(i,:) = ( D_ip2 * ( w_q(i+1,:) - w_q(i,:) ) - ...
                D_im2 * ( w_q(i,:) - w_q(i-1,:) ) )/L_node(i)^2;
            
			% d/dx( D * (dg/dx) )
            dDdg_dx2(i,:) = ( D_ip2 * ( g_q(i+1,:) - g_q(i,:) ) - ...
                D_im2 * ( g_q(i,:) - g_q(i-1,:) ) )/L_node(i)^2;
            
			% dr/dx (used in the diffusion term)
            dr_dx = ( r_q(i+1,:) - r_q(i-1,:) )/(2 * L_node(i));
            
			% C is the diffusion term needed later
            C(i,:) = w_q(i,:) .* D_i .* dr_dx.^2;
            
        else % boundary points at top and bottom (i = 1, i = N_full +1)
		
            if i == 1
                % bottom
				% use same expression as above, but define w and g at i-1 to be zero
				% this approximates the wall -> no diffusion can go through it
               
                w_q_bottom = zeros(size(w_q(i,:)));
                g_q_bottom = zeros(size(w_q(i,:)));
                
                D_ip2 = 0.5*( D_zz(i+1) + D_zz(i) ); 
                D_im2 = D_zz(i); % just take the value at i
                D_i = D_zz(i); % ditto
                
                dDdw_dx2(i,:) = ( D_ip2 * ( w_q(i+1,:) - w_q(i,:) ) - ...
                    D_im2 * ( w_q(i,:) - w_q_bottom ) ) /L_node(i)^2;
                
                dDdg_dx2(i,:) = ( D_ip2 * ( g_q(i+1,:) - g_q(i,:) ) - ...
                    D_im2 * ( g_q(i,:) - g_q_bottom ) )/L_node(i)^2;
                
                dr_dx = ( r_q(i+1,:) - r_q(i,:) )/( L_node(i)); % use one sided difference here (was central above)
                
                C(i,:) = w_q(i,:) .* D_i .* dr_dx.^2;
                
                % an alternative expression I tried out (one sided differences)
                %             D_i = D_zz(i);
                %             dDdw_dx2(i,:) = D_i * (w_q(i,:) - 2*w_q(i+1,:) + w_q(i+2,:) )/(L_node(i)^2);
                %             dDdg_dx2(i,:) = D_i * (g_q(i,:) - 2*g_q(i+1,:) + g_q(i+2,:) )/(L_node(i)^2);
                %             dr_dx = (r_q(i+1,:) - r_q(i,:) )/L_node(i);
                %             C(i,:) = w_q(i,:) .* D_i .* dr_dx.^2;
            else
                % top point (N_full + 1)
                if i >= 3
                    % enough points that I can use one sided differences
                    D_i = D_zz(i);
                    dDdw_dx2(i,:) = D_i * (w_q(i,:) - 2*w_q(i-1,:) + w_q(i-2,:) )/(L_node(i)^2);
                    dDdg_dx2(i,:) = D_i * (g_q(i,:) - 2*g_q(i-1,:) + g_q(i-2,:) )/(L_node(i)^2);
                    dr_dx = (r_q(i,:) - r_q(i-1,:) )/L_node(i);
                    C(i,:) = w_q(i,:) .* D_i .* dr_dx.^2;
                    
                else
                    % not enough points -> give up and set to 0. Means that there's only 2 points left in the tank -> can't really compute diffusion in that case
                    dDdw_dx2(i,:) = zeros(1,N_ab);
                    dDdg_dx2(i,:) = zeros(1,N_ab);
                    C(i,:) = zeros(1,N_ab);
                end
            end
        end
        
    else
        % N_full = 0
        % no diffusion
        dDdw_dx2(i,:) = zeros(1,N_ab);
        dDdg_dx2(i,:) = zeros(1,N_ab);
        C(i,:) = zeros(1,N_ab);
    end
    
    
	
    % 1st order derivatives for convection terms
	
	
	
    % MUSCL
	% uses minmod flux limiter
	% L refers to the left (bottom) side of a node boundary, R is the right (top)
	% see thesis for more details
    if i > 2 && i < N_full
        % interior grid points
        % have to be able to do i+2 and i-2
        % therefore requires: i >= 3 and N_full >= 4
        
        for j = 1:N_ab
            
            % vector of the abscissas and weights
            U = [w_q(:,j)'; g_q(:,j)'];
            
            % flow velocity
            u = u_vec(:,j)';
            
            % flux limiters at i
            r_U = (U(:,i) - U(:,i-1))./(U(:,i+1) - U(:,i));
            r_u = (u(i) - u(i-1))/(u(i+1) - u(i));
            
            flux_limiter_ui = max([0, min(1, r_u)]);
            flux_limiter_Ui(1) = max([0, min(1, r_U(1))]);
            flux_limiter_Ui(2) = max([0, min(1, r_U(2))]);
            
			% flux limiters at i + 1
            r_U = (U(:,i+1) - U(:,i))./(U(:,i+2) - U(:,i+1));
            r_u = (u(i+1) - u(i))/(u(i+2) - u(i+1));
            
            flux_limiter_uip1 = max([0, min(1, r_u)]);
            flux_limiter_Uip1(1) = max([0, min(1, r_U(1))]);
            flux_limiter_Uip1(2) = max([0, min(1, r_U(2))]);
            
			% flux limiters at i - 1
            r_U = (U(:,i-1) - U(:,i-2))./(U(:,i) - U(:,i-1));
            r_u = (u(i-1) - u(i-2))/(u(i) - u(i-1));
            
            flux_limiter_uim1 = max([0, min(1, r_u)]);
            flux_limiter_Uim1(1) = max([0, min(1, r_U(1))]);
            flux_limiter_Uim1(2) = max([0, min(1, r_U(2))]);
            
            % flux through top(i + 1/2)
            uL = u(i) + 0.5*flux_limiter_ui*(u(i+1) - u(i));
            uR = u(i+1) - 0.5*flux_limiter_uip1*(u(i+2) - u(i+1));
            
            UL = U(:,i) + 0.5*flux_limiter_Ui'.*(U(:,i+1) - U(:,i));
            UR = U(:,i+1) - 0.5*flux_limiter_Uip1'.*(U(:,i+2) - U(:,i+1));
            
            a = max([abs(uL), abs(uR)]);
            
            FR = uR*UR;
            FL = uL*UL;
            
            F_star_top = 0.5*( (FR + FL) - a*(UR - UL) );
            
            % flux through bottom(i - 1/2)
            uL = u(i-1) + 0.5*flux_limiter_uim1*(u(i) - u(i-1));
            uR = u(i) - 0.5*flux_limiter_ui*(u(i+1) - u(i));
            
            UL = U(:,i-1) + 0.5*flux_limiter_Uim1'.*(U(:,i) - U(:,i-1));
            UR = U(:,i) - 0.5*flux_limiter_Ui'.*(U(:,i+1) - U(:,i));
            
            a = max([abs(uL), abs(uR)]);
            
            FR = uR*UR;
            FL = uL*UL;
            
            F_star_bot = 0.5*( (FR + FL) - a*(UR - UL) );
			
			% flux gradient
            
            dF_dx = (F_star_top - F_star_bot)/L_node(i);
            
            duw_dx(i,j) = dF_dx(1);
            dug_dx(i,j) = dF_dx(2);
            
        end
        
    else % top and bottom boundary points
        
        if i == 1
            % bottom point
            % there's flux out the top and nothing out the bottom
			% use a 0 at the bottom, and just use upwinding at the top
            % (assuming no u bulk!!)
            
			% deal with w first
            flux_bot = 0;
            % flux_bot = -w_q(i,:) * u_bulk; % if there were a u_bulk
            flux_top = w_q(i,:).*0.5.*(u_vec(i+1,:) + u_vec(i,:));
            duw_dx(i,:) = (flux_top - flux_bot)/L_node(i);
            
			% now deal with g
            flux_bot = 0;
            % flux_bot = -g_q(i,:) * u_bulk; % if there were a u_bulk
            flux_top = g_q(i,:).*0.5.*(u_vec(i+1,:) + u_vec(i,:));
            dug_dx(i,:) = (flux_top - flux_bot)/L_node(i);
            
            
			% alternative expression: one sided differences
            %             duw_dx(i,:) = (uw_vec(i+1,:) - uw_vec(i,:))/L_node;
            %             dug_dx(i,:) = (ug_vec(i+1,:) - ug_vec(i,:))/L_node;
        else
            if i == N_full + 1
                % top point
				% now use upwinding for bottom flux (assumes velocity is upward)
				% and use values at i for top flux
				% alternative expression: backwards difference
                %                 duw_dx(i,:) = ( uw_vec(i,:) - uw_vec(i-1,:) )/(L_node);
                %                 dug_dx(i,:) = ( ug_vec(i,:) - ug_vec(i-1,:) )/(L_node);
                
				% w first
                flux_bot = w_q(i-1,:).*0.5.*(u_vec(i,:) + u_vec(i-1,:));
                flux_top = w_q(i,:).*u_vec(i,:);
                duw_dx(i,:) = (flux_top - flux_bot)/L_node(i);
                
				% now g
                flux_bot = g_q(i-1,:).*0.5.*(u_vec(i,:) + u_vec(i-1,:));
                flux_top = g_q(i,:).*u_vec(i,:);
                dug_dx(i,:) = (flux_top - flux_bot)/L_node(i);
				
            else 
                % i = 2, 3 when N_full >= 3
                % i = 2, when N_full = 2
                
                % 1st order upwind
                
                % %      rising: flux in is from -x
				
				% w first
                flux_bot = w_q(i-1,:).*0.5.*(u_vec(i,:) + u_vec(i-1,:));
                flux_top = w_q(i,:).*0.5.*(u_vec(i+1,:) + u_vec(i,:));
                duw_dx(i,:) = (flux_top - flux_bot)/L_node(i);
                
				% now g
                flux_bot = g_q(i-1,:).*0.5.*(u_vec(i,:) + u_vec(i-1,:));
                flux_top = g_q(i,:).*0.5.*(u_vec(i+1,:) + u_vec(i,:));
                dug_dx(i,:) = (flux_top - flux_bot)/L_node(i);
                
				% alternative expressions: backwards differences
                %             if i == 2
                %                 duw_dx(i,:) = (uw_vec(i,:) - uw_vec(i-1,:))/L_node;
                %                 dug_dx(i,:) = (ug_vec(i,:) - ug_vec(i-1,:))/L_node;
                %             else
                %                 % top 2 points
                %                 %             duw_dx(i,:) = ( uw_vec(i,:) - uw_vec(i-1,:) )/(L_node);
                %                 %             dug_dx(i,:) = ( ug_vec(i,:) - ug_vec(i-1,:) )/(L_node);
                %                 duw_dx(i,:) = zeros(size(uw_vec(i,:)));
                %                 dug_dx(i,:) = zeros(size(ug_vec(i,:)));
                %             end
            end
        end
    end
    
end


% check if it's become subcooled (usually a bad sign, so warn user)
if constants.step == 1 && (deltaT_sup < 0 && constants.t > 0.1)
    disp('subcooled')
end


% loop through each node with any liquid in it
for i = 1:N_full + 1
        
    % a little error checking
    if sum(abs(imag([r_q(i,:); w_q(i,:)]))) > 0
        fprintf('imaginary abscissas or weights. moments:')
        fprintf('%0.6g\t',mom)
        fprintf('\n')
        
    end
	
	
	
	
    
	% compute superheat at current node, taking into account the effect of gravity superheat
    
    % depth from surface
    depth = sum(node_level(i:N_full+1).*L_node(i:N_full+1));
    
	% gravity head
    dP_depth = rho_l*g*depth;
    
	% use curve fit I made of Tsat vs P to get effect of gravity head on superheat
	% Here's the fit: (T_sat = 19.6426*(P/1e3)^0.2499 + 122.3663) (P in kPa)
    dT_sat_depth = 19.6426*( ((P + dP_depth)/1e3)^0.2499 - (P/1e3)^0.2499 );
    
	% superheat at the current node
    deltaT_sup_node = deltaT_sup - dT_sat_depth;
	
	
	
	
	
    
    % for boiling/bubble growth, only compute if there's a decent superheat
    if deltaT_sup_node > 1e-4
	
        % Jakob number
        Ja_T = Cp_l * rho_l * C_dTs * deltaT_sup_node/(rho_tg_v * h_lv);
        
        % bubble radius rate of change
        
        % growth rate for a bubble at rest in an infinite fluid
        % simplified model: plesset & zwick
        % more complicated: scriven
        % rdot_rest = C_rdot * Ja_T^2 * alpha_l ./ r_q(i,:); % bubble at rest
        
        % equation 47 from scriven
        beta_47 = sqrt(0.5*Ja_T./(1 + (Cp_l - Cp_tg)/Cp_l * rho_tg_v/rho_l * Ja_T));
        
        % eq 47 + eq 49 from scriven
        beta_rdot = 0.85*(beta_47 + sqrt(12/pi)*beta_47.^2);
        
		% growth rate based on scriven's beta
        % use r_nuc (current value) as a crude approx of r when the bubble first nucleated
        rdot_rest = 2*beta_rdot.^2*alpha_l./(r_q(i,:));
                
        % bubble rising in the fluid. from legendre 1998
        rdot_rise = Ja_T * sqrt( 2 * alpha_l * (u_rise(i,:) + 1e-6)./...
            (pi * r_q(i,:) ) ); % rising in the liquid
        
		% sum growth rates for stationary and moving bubble
        rdot = rdot_rest + rdot_rise;
        
        % length of liquid node volume [m]
        L_l = node_level(i)*L_node(i);
        
        % surface area of node[m^2] (assumes cylindrical tank)
        A_l = pi * D_tank * L_l;
        
		% If at the bottom of the tank, include the area of the base (assume cylindrical tank)
		% Not sure why this is commented out!
        %         if i == 1
        %             A_l = A_l + 0.25*pi*D_tank^2;
        %         end
        
		% switch for different nucleation models
		% SJ: Shin & Jones (also includes other submodels within it)
		% AL: Alamgir & Lienhard (no submodels)
        switch constants.nuc_model
            
            case 'SJ' 
                
                
                % departure diameter [m]
                % correlation from Jensen & Memmel, 1986
                % gives values on the order of 10-20 microns
				% they gave 2 different correlations: one for times when you know wall superheat, one for when you don't
                switch constants.r_dep_expression
                    case 'without superheat'
                        
						% departure radius
                        r_dep = 0.5 * 2.97e4 * (P/P_cr)^-1.09 * ( K_b * T_cr / (P_cr * MW) )^(1/3);
                        
                    case 'with superheat'
                        
						% Jakob number of the wall - ambiguous in my case if it should be defined based on wall superheat, or liquid superheat
						% I ended up going with the liquid superheat.
                        %  Ja_wall = (T_lw(1) - (T_s + dT_sat_depth)) * rho_l * Cp_l/(rho_tg_v*h_lv);
                        Ja_wall = Ja_T;
						
						% constant. In my lab notebook I documented how a lot of people use different expressions for this
                        K_1 = (Ja_wall/Pr_l)^2*( g*rho_l*(rho_l-rho_tg_v)/mu_l^2 * (sigma/(g*(rho_l-rho_tg_v)))^(3/2) )^(-1);
                        
						% Eotvos number at departure
                        Eo_dep = ( 0.19*( 1.8 + 1e5*K_1)^(2/3) )^2;
						
                        % radius at departure
                        r_dep = 0.5 * (Eo_dep*sigma /(g*(rho_l - rho_tg_v)) )^0.5;
                        
                    otherwise
                        error('invalid constants.r_dep_expression. try ''without superheat'' or ''with superheat''')
                end
                
                
                
                % radius of new bubbles
                % common expression
                %         r_nuc = C_r_nuc * 2*sigma*T_s/(rho_tg_v * h_lv * C_dTs * deltaT_sup_node);
                % more advanced one. not sure where I got it from, but it seems to
                % be equal to about 2x the above (2.1 or 2.2 generally)
                % it's listed in hibiki & ishii 2003 paper they also show how it
                % reduces to the common expression
                r_nuc = (2*sigma*(1 + rho_tg_sat/rho_l)/P)...
                    /( exp( h_lv * (deltaT_sup_node)/(Ru/MW * T_l*T_s)) - 1);
                
                
                % hysteresis!
                % see Qi & Klausner, 2005
                % basically, if we're depressurizing, the meniscus is
                % concave and you need 2x the superheat to activate it
                if constants.include_hysteresis
                    if deltaT_sup_node < constants.deltaT_sup_max
                        % we're past max superheat
                        % so now the menisucs is convex
                        r_nuc_max = (4*sigma*(1 + rho_tg_sat/rho_l)/P)...
                            /( exp( h_lv * (constants.deltaT_sup_max)...
                            /(Ru/MW * T_l*T_s)) - 1);
                        r_nuc = max([ r_nuc, r_nuc_max ]);
                    else
                        % we haven't reached max superheat yet
                        r_nuc = 2*r_nuc;
                    end
                end
                
                % nucleation site density
                switch constants.nuc_density_expression
                    
                    case 'shin and jones'
                        % nucleation density (shin and jones 1993) and (blinkov,
                        % jones, nigmatulin, 1993)
                        
                        % non-dimensional cavity size (ie bubble nucleation size)
                        r_c_star = r_nuc / r_dep;
                        
                        % non-dimensional nucleation rate
                        N_ns_star = 1e-7 * r_c_star^-4;
                        
                        % nucleation site density [1/m^2]
                        nuc_density = N_ns_star * ( 0.5 / r_dep )^2; % original expression
                        
                        % other ones I've tried
                        %                 nuc_density = N_ns_star * ( 0.5 / r_nuc )^2;
                        %                 nuc_density = 1e-4 * ( 0.5 / r_nuc )^2;
                        
                    case 'hibiki and ishii'
                        % fancy expression for nucleation density
                        % (hibiki & ishii, 2003)
                        
                        % originally thought:
                        % gives much larger values than the shin/jones one
                        % roughly 10^10 larger, with actual values around 10^8
                        
                        % then found I left out the ^-4 on the N_ns_star relation
                        % givens much smaller values, by about 10^-6
                        % actual values around 10^5
                        
                        % depends on the contact angle - I found a paper where
                        % they measured contact angle of CO2 on SS316 and made a
                        % curve fit from their results
                        N_nbar = 4.72e5; % constant
                        mu_HI = 0.722; % constant
                        lambda_prime = 2.5e-6; % constant
                        % advancing and receding contact angles (in degrees)
                        ACA = -0.003417*(T_s - 273.15)^2 - 0.2873*(T_s - 273.15) + 29.83;
                        RCA = -0.004171*(T_s - 273.15)^2 - 0.3386*(T_s - 273.15) + 16.38;
                        theta_HI = deg2rad(mean([ACA RCA])); % take contact angle as mean of advancing and receding
                        
                        R_c = r_nuc; % just matching their notation
                        
                        rho_plus = log10(delta_rho/rho_tg_sat); % non-dimensional delta rho
                        
                        f_rho_plus = -0.01064 + 0.48246*rho_plus - 0.22712*rho_plus^2 + 0.05468*rho_plus^3; % empirical fit
                        
                        nuc_density = N_nbar *( 1 - exp( - theta_HI^2/(8*mu_HI^2) ) )...
                            *(exp( f_rho_plus*lambda_prime/R_c) - 1); % nucleation rate (#/m^2)
                    otherwise
                        error('invalid constants.nuc_density_expression. try ''shin and jones'' or ''hibiki and ishii''')
                end
                
                % nucleation frequency
                switch constants.nuc_frequency_expression
                    case 'shin and jones'
                        
                        % nucleation frequency (shin and jones also)
                        
                        % nucleation frequency [Hz]
                        nuc_freq = 1e4 * C_dTs * deltaT_sup_node^n_nuc_freq;
                        
                    case 'saddy and jameson'
					% it's really just proportional to r^2, but I have other constants that can take care of that
                        nuc_freq = r_dep^2;
                    otherwise
                        error('invalid constants.nuc_frequency expression. try ''shin and jones'' or ''saddy and jameson''')
                end
                
                % nucleation rate [Hz]
                nuc_rate = C_nuc_rate * nuc_density * nuc_freq * A_l;
				
            case 'AL'
                
                % alamgir and lienhard, 1981:
				% more similar to a homogeneous nucleation model
                
                B = K_b * T_l / h_planck;
                
                v_g = 1/rho_tg; % specific volume of gas
                v_f = 1/rho_l; % specific volume of liquid
                
                J = ( N_A/( 1/rho_l * MW) )^(2/3) * B * ...
                    exp( -16*pi*sigma^3*phi / ( 3*K_b*T_s*(1 - v_f/v_g)^2 * (P_sat - P)^2 ) ); % nucleation rate per area
                
                nuc_rate = J * A_l;
                
        end
        
        % homogeneous nucleation rate
        
        rstar = 2*sigma/(P_sat - P);
        
        N_l = rho_l / MW * 1e3*N_A;
        
        m_M = 1e3 * MW / N_A;
        
        J_hom = N_l * rho_l/rho_tg * sqrt(2*sigma/(pi * m_M)) ...
            * exp( - 4*pi*sigma*rstar^2/(3*K_b*T_l));
        
        nuc_rate = nuc_rate + J_hom * V_node(i);
        
        if constants.step == 1 && (J_hom * V_node(i)/nuc_rate > 1e-6)
            disp('we''ve got some homogeneous nucleation happening')
        end
        
        
    else
        
        % no superheat -> no bubble growth or nucleation
        r_nuc = 0;
        rdot = 0;
        nuc_rate = 0;
        
    end
    
    % rdot term based on rho_dot
    rdot_rhodot = 0;
    if constants.include_rdot_rhodot
        if ~isnan(guesses.rhodot_tg_sat)
            rdot_rhodot = -guesses.rhodot_tg_sat/rho_tg_sat * r_q(i,:)/3;
        end
    end
    
    % generat scaled abscissas for more accurate calculations
    % later on, _s will signify scaled (I think)
    r_m = max(abs(r_q(i,:)));    % value to scale by
    
    r_s = r_q(i,:)/r_m;  % scaled abscissas
    
    % set the first term to 0 because growth has no effect on 0th moment
    % (total # of bubbles) see page 71 of notebook #5 for details
    growth_int_s = zeros(1,2*N_ab);
    
    % calculate the scaled growth integral term
    for k = 2:N_ab*2
        growth_int_s(k) = ((k-1)/p)*sum( (1/r_m) * r_s.^((k-1)/p - 1) .* w_q(i,:) .* (rdot + rdot_rhodot));
    end
    
    % nucleation rate per volume (specific nucleation rate)
    
    % think this maybe shouldn't have the node level term in it, and
    % maybe not the V node either -> maybe V_l_star?
    % changed my mind -> nuc_rate is based on V_node
    spec_nuc_rate = nuc_rate / ( (node_level(i) + 1e-3) * V_node(i) );
	% add the 1e-3 just in case node_level is very low and creates a huge value when inverted
    

	% some printing to screen for debugging
    %     if (i == 1 && constants.step == 1)
    %         fprintf('spec nuc rate = %0.4g\n growth_int = %0.4g\n',spec_nuc_rate,growth_int_s(V_moment_index)*r_m^3)
    %     end
    
    
    % preallocate
    birth_int_s = zeros(1,N_ab*2);
    
	% if we have a decent superheat, calculate birth integral term
    if deltaT_sup_node > 1e-4
        
        for k = 1:N_ab*2
		% the commented code here is for different birth distributions
		
            % if nucleation happens only at r_nuc (delta fn distribution)
            birth_int_s_delta = (r_nuc/r_m).^((k-1)/p) * spec_nuc_rate;
            
            %      % exponential distribution
            %      r_a = 10*r_nuc;
            %
            %      birth_int_s_exp = (r_a/r_m)^((k-1)/p) * spec_nuc_rate * exp( r_nuc/r_a ) ...
            %        * gamma(1+((k-1)/p)) * gammainc(r_nuc/r_a, 1+((k-1)/p), 'upper');
            %                     %
            %         uniform distribution
            %         dr = 100*r_nuc;
            %         birth_int_s_uni = (1/r_m)^((k-1)/p) * spec_nuc_rate/dr * 1/( (k-1)/p + 1)*...
            %                     ( (r_nuc + dr)^( (k-1)/p + 1) - r_nuc^( (k-1)/p + 1) );
            
            
            birth_int_s(k) = birth_int_s_delta;
            
            % now if it's spread out a bit (gaussian distribution)
            %         s_nuc = r_nuc;
            % %         birth_int_s(k) = spec_nuc_rate * r_m*integral(@(rs) rs.^k .* ...
            % %             (1/(s_nuc*sqrt(2*pi))).*exp( - ((r_m*rs) - r_nuc).^2./(2*s_nuc^2)), 0, 25*r_nuc/r_m);
            %
            %         rs = linspace(0,r_nuc/r_m + 5*s_nuc/r_m,1000);
            %         birth_int_s(k) = spec_nuc_rate * r_m * trapz(rs, rs.^((k-1)/p) .* ...
            %             (1/(s_nuc*sqrt(2*pi))).*exp( - ((r_m*rs) - r_nuc).^2./(2*s_nuc^2)));
        end
        
    end
    
    % birth and death due to coalescence
    
    coal_birth_s = zeros(2*N_ab,1);
    coal_death_s = zeros(2*N_ab,1);
    
    % only bother if I've turned it on
    if C_coalescence(1) > 0
        
        % most of this is from prince and blanch (1990)
        % buoyancy and turbulence driven
        
        nu_t = 0.0536*D_tank^1.77/rho_l; % turbulent viscosity
        U_max = ( (1 - 0.75*V_bubi(i))/(1 - V_bubi(i)) )*V_bubi(i) * D_tank^2/(48*nu_t); % max of mean velocity expected in centerline
        mean_shear = 0.53*U_max/(0.5*D_tank); % mean shear from mean velocity
        
        P2 = P; % used for calculating pressure difference
        
        liquid_height = V_l_star/(pi*D_tank^2/4); % height of liquid level
        
        P1 = P2 + rho_l*g*liquid_height; % gravity head (there's something fishy here -> should probably be rho of mixture)
        
        Q = V_tank/8; % volumetric gas flow rate. just assumed a constant value here
        
        turb_diss = Q*g * P2*log(P1/P2) / ( pi * (0.5*D_tank)^2 * (P1 - P2) ); % turbulent dissipation rate
        
		% loop through and calculate coalescence rates for each combination of bubble size
        for k = 1:N_ab*2
            for l = 1:N_ab
                
                
                % radius and diameter of bubble i and bubble j
                % that was the original indices, had to change when I
                % went to a 1D model. now it's really l and j
                % j is the vector here, summed over j at the end
                rbi = abs(r_q(i,l));
                rbj = abs(r_q(i,:));
                dbi = 2*rbi;
                dbj = 2*rbj;
                
                % ----- theta for laminar shear -----
                qLS = C_coalescence(2) * 4/3*(rbi + rbj).^3*mean_shear;
                
                % rise velocity for bubble i and j
                %                 u_ri = sqrt( (2.14*sigma/(rho_l*dbi)) + 0.505*g*dbi);
                %                 u_rj = sqrt( (2.14*sigma./(rho_l*dbj)) + 0.505*g*dbj);
                u_ri = u_rise(i,l);
                u_rj = u_rise(i,:);
                
                % ----- collision area -----
                Sij = pi/4*(rbi + rbj).^2;
                
                % ----- theta for buoyancy -----
                qB = Sij.*abs(u_ri - u_rj);
                
                % ----- theta for turbulence -----
                qT = C_coalescence(3) * 0.089*pi*(dbi + dbj).^2 .* turb_diss.^(1/3) .* sqrt(dbi.^(2/3) + dbj.^(2/3));
                
                % ----- equivalent bubble radius -----
                rb_eq = ( (1./rbi + 1./rbj)/2 ).^-1;
                
                % ----- contanct time -----
                
                % prince + blanch, 1990, in turn from Levich, 1962
                % note I added the 0.1 factor in front. P+B say this is
                % really just an order of magnitude (if that) estimate
                t_cont = 0.1*rb_eq.^(2/3) ./ turb_diss.^(1/3);
                
                % kamp & chesters, 2001
                %                 rho_c = rho_l;
                %                 C_vm = 0.8;
                %                 t_cont = sqrt( rho_c*C_vm/(3*sigma) * ( 2*dbi.*dbj./(dbi + dbj)).^3 );
                
                % film initial and final thicknesses
                film_i = 1e-4; % wild guess
                film_f = (C_hamaker * rb_eq./(8*pi*sigma)).^(1/3);
                
                % time required for coalescence
                t_coal = sqrt( rb_eq.^3 * rho_l/(16 * sigma) ) .* log( film_i ./ film_f);
                
                % coalescence kernel
                beta = C_coalescence(1)*(qT + qB + qLS).*exp( - t_coal ./ t_cont);
                
                % add to net coalescence rates
                coal_birth_s(k) = coal_birth_s(k) + sum( 0.5 * w_q(i,l) * w_q(i,:) .* ( abs(r_s(l))^3 + abs(r_s).^3 ).^((k-1)/(3*p)) .* beta );
                coal_death_s(k) = coal_death_s(k) + sum( r_s(l)^((k-1)/p) * w_q(i,l) * w_q(i,:) .* beta );
                
            end
        end
        
        
    end
    
	% check to make sure that coalescence is conserving mass
    if abs(coal_birth_s(V_moment_index) - coal_death_s(V_moment_index))/coal_birth_s(V_moment_index) > 1e-6
        disp('coalescence isn''t conserving mass')
    end
    
	% rate of change of the moments (scaled)
    dmom_dt_s = birth_int_s(:) + growth_int_s(:) + coal_birth_s(:) - coal_death_s(:);
	
	% switch to notation that I found in marchisio paper
    beta_q = dmom_dt_s;
    
	% initialize several things
    not_converging = 1;
    linear_eqn_counter = 0;
    r_sp = r_s;
    
    C_star = C(i,:)'/r_m^2; % scale
    
	% iterate to converge on the Ax = b solution
	% (here x = alpha, b = d)
	% when it has problems, I will jiggle the abscissas a little bit
    while not_converging
        
        % pre allocate the parts of A
        A1 = zeros(2*N_ab, N_ab);
        A2 = A1;
        A3 = A1;
		
		% define A
        for j = 0:(2*N_ab - 1)
            
            if j == 0
                A1(1,:) = ones(1,N_ab);
                A2(1,:) = zeros(1,N_ab);
                A3(1,:) = zeros(1,N_ab);
            elseif j == 1
                A1(2,:) = (p-1)/p * r_sp.^(1/p);
                A2(2,:) = (1/p) * r_sp.^(1/p - 1);
                A3(2,:) = zeros(1,N_ab);
            elseif j == 2
                A1(j+1,:) = (1 - j/p) * r_sp.^(j/p);
                A2(j+1,:) = (j/p) * r_sp.^(j/p - 1);
                A3(j+1,:) = 2*ones(1,N_ab);
            else
                A1(j+1,:) = (1 - j/p) * r_sp.^(j/p);
                A2(j+1,:) = (j/p) * r_sp.^(j/p - 1);
                A3(j+1,:) = j*(j-1) * r_sp.^(j - 2);
            end
        end
        
        A = [A1 A2];
        
        d = A3*C_star + beta_q;
        
		% solve for alpha
        [alpha_q, error_flag] = linear_equation_solver(A,d);
        
		% calculate the relative spacing between the abscissas
		% if it's too small, can cause problems (this is only used for printing to screen)
        dr = diff(r_s);
        for k = 1:length(dr)
            rbar = mean([r_s(k) r_s(k+1)]);
        end
		
		% if there was a problem in the linear equation solver, try jiggling the abscissas a little bit
        if error_flag == 1
		
            r_sp = r_s + (rand(size(r_s)) - 0.5).*abs(r_s)*1e-6;
			
            linear_eqn_counter = linear_eqn_counter + 1;
			
			% if we've jiggled a bunch of times and it's not working, just use "\" and give up
            if linear_eqn_counter >= 25
                fprintf('linear equation solver not converging, rcond = %0.4g, min dr/r = %0.4g\n',rcond(A), min(abs(dr./rbar)) )
                alpha_q = A\d;
                not_converging = 0;
            end
            
        else
			% no errors, we've converged
            not_converging = 0;
            
        end
        
    end
	
	% extract a and b* from alpha
    a_q = alpha_q(1:N_ab);
    b_q_s = alpha_q(N_ab+1:end);
	
    b_q = b_q_s * r_m; % go from scaled b (ie b*) back to the real thing
    
    % have to put this in the right place
    % for i = 1, it should be 1, 2, ... N_ab
    % for i = 2, N_ab+1, N+2, ... N_ab+N_ab
    % for i = i, it's (i-1)*N_ab, ... i*N_ab;%
    ind_node = (i-1)*N_ab + [1:N_ab];
    
    
	% rate of change of the weights and weighted abscissas
	
    % only include birth/death/growth (ie 0D)
    %         dw_dt(ind_node) = a_q;
    %         dg_dt(ind_node) = b_q;
    
    % include flux of bubbles in physical space (including diffusion)
    dw_dt(ind_node) = a_q - duw_dx(i,:)' + dDdw_dx2(i,:)';
    dg_dt(ind_node) = b_q - dug_dx(i,:)' + dDdg_dx2(i,:)';
    
	% birth and growth terms
	% extract from the moments and undo the scaling
    % these have units of (m^3/(m^3*s)) ie (volume/time)/(volume of mix)
    birth_term(i) = birth_int_s(V_moment_index)*r_m^(3);
    growth_term = growth_int_s(V_moment_index)*r_m^(3);
    
	% bubbles leaving the free surface and entering the ullage
	% only if we're at the node with the liquid surface in it
	% note units of the death term: m^3/(m^2 * s)
    if i == N_full + 1
    
        %  if we're looking at the bottom node, just take its value
        if i == 1
            death_term = 4/3 * pi * sum(r_q(i,:).^(3) .* w_q(i,:) .* (u_rise(i,:) - u_LL) );
        else
		
            % above the bottom node, linearly interpolate/extrapolate to
            % get the value wherever the free surface is
			% the commented code is from other schemes I tried out
            death_term_i = 4/3 * pi * sum(r_q(i,:).^(3) .* w_q(i,:) .* (u_rise(i,:) - u_LL) );
            death_term_im1 = 4/3 * pi * sum(r_q(i-1,:).^(3) .* w_q(i-1,:) .* (u_rise(i-1,:) - u_LL) );
            % %             if node_level(i) < 0.5 % only interpolate if the liquid level is below 50% in a node
            % %                 death_term = (0.5 + node_level(i))*death_term_i + (0.5 - node_level(i))*death_term_im1;
            % %             else
            %             death_term_slope = (death_term_i - death_term_im1)/1;
            %             death_term = death_term_i + death_term_slope*( node_level(i) - 0.5 );
            death_term_slope = (death_term_i - death_term_im1)/(L_node(i)/2 + L_node(i-1)/2);
            death_term = death_term_im1 + (L_node(i-1)/2 + node_level(i)*L_node(i)) * death_term_slope;
            % %             end
            % %             else % use a parabola for fitting instead of a line
            % %                death_term_i = 4/3 * pi * sum(r_q(i,:).^(3) .* w_q(i,:) .* (u_rise(i,:) - u_bulk) );
            % %                death_term_im1 = 4/3 * pi * sum(r_q(i-1,:).^(3) .* w_q(i-1,:) .* (u_rise(i-1,:) - u_bulk) );
            % %                death_term_im2 = 4/3 * pi * sum(r_q(i-2,:).^(3) .* w_q(i-2,:) .* (u_rise(i-2,:) - u_bulk) );
            % %                death_term = interp1( [-2 -1 0], [death_term_im2, death_term_im1, death_term_i], (node_level(i)-0.5),'nearest','extrap');
            % %             end
        end
    else
	% if not at the top, there's no bubble death
        death_term = 0;
    end
    
	% at the bottom of the tank, have to worry about bubbles being convected out of the tank
	% I zeroed this out because I decided this should only be used if the bubbles were leaving at a different rate than the liquid (ie composition changing)
    if i == 1
        %         death_term_injector = 4/3 * pi * sum(r_q(i,:).^(3) .* w_q(i,:) .* (-u_rise(i,:) + u_bulk) );
        %         if death_term_injector < 0
        death_term_injector = 0;
        %         end
    end
    
    % mdot into bubbles from liquid
    mdot_bub_l(i) = node_level(i) * V_node(i) * 4/3 * pi * rho_tg_sat * (birth_term(i) + growth_term - (death_term/(4/3*pi)));
    
    % mdot into bubbles from tg (really from bubbles into tg, but have to
    % keep sign convention for mdot)
    mdot_bub_tg(i) = - pi/4*D_tank^2 * rho_tg_sat * death_term;
    
    % mdot into the bubbles
    mdot_bub(i) = mdot_bub_l(i) + mdot_bub_tg(i);
    
	% error checking
    if isinf(mdot_bub)
        disp('inf problem')
    end
    
end

% in the empty nodes, there should be no mass going in/out of bubbles
% (because there aren't any bubbles there)
mdot_bub_tg(N_full+2:N_nodes) = 0;
mdot_bub_l(N_full+2:N_nodes) = 0;

% mass flow rate of bubbles out the injector
mdot_bub_injector = pi/4*D_tank^2 * rho_tg_sat * death_term_injector;

% need to deal with dw_dt and dg_dt for the nodes that have left the liquid
% set them equal to the top node in case the fill level rises
% (don't want a discontinuity in that case)
for i = N_full+2:N_nodes
    
    ind_node = (i-1)*N_ab + [1:N_ab];
    
    ind_Nfp1 = ((N_full + 1) -1)*N_ab + [1:N_ab]; % index of N_full + 1
    
    dw_dt(ind_node) = dw_dt(ind_Nfp1);
    dg_dt(ind_node) = dg_dt(ind_Nfp1);
    
end

% net rate of change of gas mass (mdot_bub_tg is negative)
mdot_tg = sum( -mdot_bub_tg );

% net rate of change of liquid mass
mdot_l = - sum(mdot_bub_l) - mdot_out_liq;

% HT from wall to liquid (natural convection correlation)
% Qdot_lw = Qdot('lw',T_l,T_lw(1),rho_l,m_l,D_tank);

% boiling heat transfer model from Gorenflo and Kotthoff, 2005
Pr = P/P_cr; % reduced pressure
h_20_lw = 1e3 * exp( 0.3092*log(Pr)^3 + 1.649*log(Pr)^2 + 3.641*log(Pr) + 5.272); % heat transfer coeff. at 20 kW/m^2
n_lw = 1 - 0.3*Pr^0.3; % exponent used

% reference heat flux (20 kW/m^2)
q_20 = 20e3;

% if the wall is above saturation temperature, calculate heat flux based on boiling correlation
if T_lw(1) > T_s
    q_lw = ( ( (T_lw(1) - T_s)*h_20_lw )^(1/n_lw) / q_20 ) ^(n_lw/(1 - n_lw));
else
    q_lw = 0;
end

% wetted tank wall area
A_l = 4*V_l/D_tank + pi/4*D_tank^2;

% HT from liquid wall into liquid
Qdot_lw = C_qdot_lw * q_lw*A_l;

% net HT into liquid
Qdot_l = Qdot_lw;

% HT into gas from wall
Qdot_gw = Qdot('gw',T_tg,T_gw(1),rho_tg,m_tg,D_tank);

% net HT into gas
Qdot_tg = Qdot_gw;

% rate of change of total bubble volume
% I was never sure if this should include a term to account for the changing density
Vdot_bub = (sum( mdot_bub ) - mdot_out_vap)/ rho_tg_sat;

% rate of change of internal energy for liquid and ullage
% this isn't actually Udot, but Udot without the P*Vdot term (hence the i)
Udot_li = - mdot_out_liq*h_l - sum( mdot_bub_l )*((h_tg_sat - h_l) + h_l) + Qdot_l - P*Vdot_bub;

Udot_tgi = Qdot_tg - sum( mdot_bub_tg )*(h_tg_sat) ;

% partial derivative of specific internal energy w.r.t. density for liquid
% (used to need this)
% du_drho_tg = dh_drho_tg + P/rho_tg^2  - 1/rho_tg * dP_drho_tg;

% partial derivative of specific internal energy w.r.t. density for liquid
du_drho_l = dh_drho_l + P/rho_l^2  - 1/rho_l * dP_drho_l;

% compute the rate of change of the liquid volume
Vdot_l = solve_for_Vdot(Udot_tgi, mdot_tg, m_tg, ...
    Udot_li, u_l, mdot_l, m_l, du_drho_l, Cv_l, dP_drho_l, V_l, ...
    dP_dT_l, V_tg, P, drho_dx_P_tg, drho_dP_x_tg, u_tg_v_sat, ...
    u_tg_l_sat, x_tg, ...
    du_dT_sat_tg_v, du_dT_sat_tg_l, dP_dT_tg_sat, ...
    rho_tg_v, u_tg, guesses, Vdot_bub);


% rate of change of ullage volume
Vdot_tg = - Vdot_l - Vdot_bub;

% rate of change of ullage internal energy
Udot_tg = Udot_tgi - P*Vdot_tg;

% rate of change of liquid internal energy
Udot_l = Udot_li - P*Vdot_l;

% rate of change of liquid density
rhodot_l = mdot_l/V_l - m_l/V_l^2 * Vdot_l;

% rate of change of temperature of the liquid
Tdot_l = ( ( Udot_l - u_l*mdot_l )/m_l - du_drho_l*rhodot_l )/Cv_l;

% mass of wall exposed to liquid
m_lw = tank_wall_mass(V_l,D_tank,rho_w,t_w);

% mass of wall exposed to gas
m_gw = tank_wall_mass(V_tg,D_tank,rho_w,t_w);

% HT from air to gas wall
Qdot_agw = Qdot('agw',T_air,T_gw(end),rho_tg,m_tg,D_tank);

% HT from air to liquid wall
Qdot_alw = Qdot('alw',T_air,T_lw(end),rho_l,m_l,D_tank);

% conduction from liquid wall to gas wall (this isn't actually used right now!)
L_tank = 4*V_tank/(pi*D_tank^2);
L_wc = L_tank/2; % length scale for wall conduction 
Qdot_wc = k_w*(T_lw - T_gw)*pi*D_tank*t_w/L_wc;

% rate of change of mass of gass wall
mdot_gw = 4*Vdot_tg*t_w*rho_w/D_tank;

% wall conduction

% first some older simple expressions for lumped model of wall
% rate of change of temperature of gas wall
% Tdot_gw = (Qdot_agw - Qdot_gw + Qdot_wc + cv_w*mdot_gw*(T_lw - T_gw))/(m_gw*cv_w);

% rate of change of temperature of liquid wall
% Tdot_lw = (Qdot_alw - Qdot_lw - Qdot_wc)/(m_lw*cv_w);

% or instead can solve transient conduction equation in wall
Tdot_lw = wall_conduction(T_lw, Qdot_lw, Qdot_alw, constants);
Tdot_gw = wall_conduction(T_gw, Qdot_gw, Qdot_agw, constants);

alpha_w = k_w/(rho_w * cv_w); % thermal diffusivity of wall

% compute fourier number for wall conduction (see if it's going unstable!)
fourier_number = alpha_w*constants.h/(t_w/constants.N_rw)^2;
% 
% if fourier_number < 0.3
%     disp('fourier number is less than 0.3')
% end

% output dy/dt
dy = [mdot_tg;
    Udot_tg;
    Tdot_gw;
    mdot_l;
    Tdot_l;
    Tdot_lw;
    dw_dt';
    dg_dt'];


% compute simple prediction of y at new time step
% if it looks like it's going negative, warn the user
% (don't want negatives in any part of y)
y_predictions = y + dy * constants.h;

if min(y_predictions) < 0
    disp('predicting y < 0')
end

if ~sum(isreal(dy(:))) || sum(isnan(dy(:)))
    disp('complex or nans in derivatives')
    error_flag = 1;
end

% probably need to check to make sure it's not going negative in nodes
% above N_full+1
% correct for negative stuff
% only if we found negative entries and dy for those is < 0
% (otherwise if dy > 0, it was correcting itself)
if ~isempty(ind_negative) && sum( dy(ind_negative) < 0 ) > 0
    % we have negative dydt where y < 0!
    for i = 1:length(ind_negative)
        if dy(ind_negative(i)) < 0
            dy(ind_negative(i)) = 0;
        end
    end
    
    error_flag = 1;
    
    disp('negative stuff in diffeqns')
end
%
% % check mass conservation
% % mass out through injector
%
% % mass of ullage
% % mass of liquid
% % mass transfer between them
% % mass of bubbles
%
% % net rate of change of gas mass (mdot_bub_tg is negative)
% mdot_tg = sum( -mdot_bub_tg );
%
% % net rate of change of liquid mass
% mdot_l = - sum(mdot_bub_l) - mdot_out_liq;
if V_bubi_error_flag ==1 && error_flag == 0
    error_flag = 1;
end

% if only 1 output, it's just dy
% if 2, also pass lots of variables
if nargout == 1
    
    varargout{1} = dy;
    
else
    varargout{1} = dy;
    
    m_bub = V_bub*rho_tg_sat;
    U_bub = m_bub*u_tg_v_sat;
    n_bubi = sum_over_nodes(mom(:,1), node_level, V_node)/V_l_star;
    fill_level = V_l_star/V_tank;
    gas_holdup = V_bub/V_l_star;
    U_liq = u_l*m_l;
    
	% store lots of variables in the debug_data structure
    debug_data.m_tg = m_tg;
    debug_data.U_tg = U_tg;
    debug_data.m_l = m_l;
    debug_data.T_l = T_l;
    debug_data.g_q = g_q;
    debug_data.r_q = r_q;
    debug_data.mom = mom;
    debug_data.w_q = w_q;
    debug_data.V_bubi = V_bubi;
    debug_data.T_s = T_s;
    debug_data.P = P;
    debug_data.deltaT_sup = deltaT_sup;
    debug_data.rho_l = rho_l;
    debug_data.rho_tg = rho_tg;
    debug_data.rho_tg_sat = rho_tg_sat;
    debug_data.x_tg = x_tg;
    debug_data.V_l = V_l;
    debug_data.V_tg = V_tg;
    debug_data.V_bub = V_bub;
    debug_data.m_bub = m_bub;
    debug_data.U_bub = U_bub;
    debug_data.V_l_star = V_l_star;
    debug_data.n_bubi = n_bubi;
    debug_data.fill_level = fill_level;
    debug_data.Qdot_lw = Qdot_lw;
    debug_data.gas_holdup = gas_holdup;
    debug_data.N_full = N_full;
    debug_data.U_liq = U_liq;
    debug_data.mdot_out_liq = mdot_out_liq;
    debug_data.mdot_out_vap = mdot_out_vap;
    debug_data.h_l = h_l;
    debug_data.h_tg_sat = h_tg_sat;
    debug_data.gas_holdup_injector = V_bubi(1);
    debug_data.Vdot_bub = Vdot_bub;
    debug_data.T_lw = T_lw;
    debug_data.T_gw = T_gw;
    debug_data.u_rise = u_rise;
    debug_data.N_bubi = N_bubi;
    debug_data.N_full = N_full;
    debug_data.mdot_tg = mdot_tg;
    debug_data.mdot_bub_l = mdot_bub_l;
    debug_data.mdot_bub_tg = mdot_bub_tg;
    debug_data.mdot_l = mdot_l;
    debug_data.Cp_l = Cp_l;
    debug_data.node_level = node_level;
    
    debug_data.diff_eqns_error_flag = error_flag;
    
    varargout{2} = debug_data;
end



