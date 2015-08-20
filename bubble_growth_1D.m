function varargout = bubble_growth_1D(varargin)
% test for github
tic
% close all
warning('off','MATLAB:nearlySingularMatrix')
warning('off','MATLAB:Axes:NegativeDataInLogAxis')
warning('off','MATLAB:interp2:NaNstrip')
%% define input parameters

switch computer
    case {'GLNXA64'}
        % on server
        addpath('~/refprop/refprop');
        t_save = 60*60; % save interval in seconds
        plot_stuff = 0;
        save_stuff = 1;
        save_periodically = 0;
        save_parameters_only = 1;
        time_out = 1;
        max_comp_time = 60*60;
        results_save_dir = pwd;
        
        
    case {'MACI64','PCWIN64'}
        % on laptop
        t_save = 3*60; % save interval in seconds
        plot_stuff = 1;
        save_stuff = 1;
        save_periodically = 0;
        save_parameters_only = 0;
        plot_periodically = 0;
        t_plot = 5;
        time_out = 1;
        max_comp_time = 6*60*60;
        if strcmp(computer,'PCWIN64')
            results_save_dir = 'C:\Users\jonah\compiled research\model_results';
        else
            results_save_dir = '/Users/jez/School/stanford/compiled research/tank modeling/model_results_data_files';
        end
        
end
current_dir = pwd;
% constants.property_source = 'PDT';
constants.fsolve_options = optimset('TolX',1e-12,'Display','off');
constants.property_source = 'refprop';
constants.nuc_model = 'SJ';
constants.ADQMOM_p = 1;
p = constants.ADQMOM_p;
fluid = 'CO2';
constants.fluid = fluid;

ADQMOM = 'off';

if nargin == 0

N_nodes = 50;
N_mom = 4;
rel_tol = 1e-3;     % [] max relative error allowed in adaptive scheme
constants.C_qdot_lw = 2e-4;
constants.C_coalescence = [0 5e-6 5e-2]; % collision efficiency, laminar shear, turbulence
constants.C_nuc_rate = 6e4;

else
    inputs = varargin{1};
    N_nodes = inputs.N_nodes;
    N_mom = inputs.N_mom;
    rel_tol = inputs.rel_tol;
    
    constants.C_qdot_lw = inputs.C_qdot_lw;
    constants.C_coalescence = inputs.C_coalescence;
    constants.C_nuc_rate = inputs.C_nuc_rate;

end

constants.N_nodes = N_nodes;

IC_moments = gamma_dist_moments( 5e-10, 1e-3, N_mom, p);

alpha_ic = 1e-12;

N_ab = N_mom/2;
constants.N_ab = N_ab; % number of abscissas
% (2*N = number of moments, going from 0 to 2N-1)

constants.phi = 0.1;

clock_save = clock;
clock_start = clock;
clock_plot = clock;

% close all
N_rw = 20;
constants.N_rw = N_rw;

constants.C_rdot = 1;%(3/(2*pi));%2.5*pi;
constants.C_u_rise = 1e0;
constants.n_nuc_freq = 3;
constants.C_r_nuc = 1;
constants.C_dTs = 1;
constants.C_x_inj = 1;
if constants.C_coalescence(1) == 0
constants.coalescence_switch = 'off';
else
    constants.coalescence_switch = 'on';
end
save_filename = 'bubble_sim_data';
d_inj = 0.0708 * 0.0254;
% d_inj = 0.056 * 0.0254;
A_inj = pi/4 * (d_inj^2);
%     A_inj = 1.8e-8;
specified_case = 12;

constants.f_feedline = 0.0;
constants.L_feedline = 5 * 0.0254;
constants.D_feedline = (.375 - 2*0.049) * 0.0254;

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
    D = sqrt(4/pi*V_tank/L_tank);
    % [m] tank diameter
    k_w = 0.195;          % [W/m.K] thermal conductivity of wall
else
    
    [Ti, fill_level, V_tank, L_tank, ...
        Cd, Po, T_air, rho_w, cv_w, t_w, D, k_w] = initial_conditions(specified_case);
    
end

% divide the whole tank into nodes
% L_node = L_tank*(1 + 20*alpha_ic)*fill_level/(N_nodes - 0.5);
L_node = L_tank/(N_nodes-1); % use N -1 if there are points on the boundaries
V_node = pi*0.25*D^2*L_node;

constants.L_node = L_node;
constants.V_node = V_node;

% initialize program parameters
h = 1e-12;           % [s] initial time step
running = 1;        % [] switch, 1 = program running, 0 = program stopped
abs_tol = 1e9;     % [] max absolute error allowed in adaptive scheme
min_error = 1e-4;   % [] min error (relative to error_tol) before step size is increased
h_max = 1;       % [s] max allowable time step
h_min = 1e-16;      % [s] min allowable time step
t_end = 1e3;         % [s] end time (if LRO doesn't happen first)
LRO_tol = 2e-1;     % [s] tolerance for resolving the LRO point
dT_sup_tol = 1e-14;

ode_solver = 'DP'; % [] options:
% 'RKF' for runge-kutta-fehlberg
% 'euler' for 1st order euler
% 'RK4' for 4th order runge-kutta
% 'CK' for cash-karp
% 'DP' for dormand-prince


% 6 variables for liquid and vapor
% then N_mom for each node (N_mom/2 abscissas, N_mom/2 weights)
N_dim = 4 + 2*N_rw + N_mom*N_nodes;

% fsolve_options = optimset('display','off');

hesson_fit = load('hesson_fit');
constants.hesson_fit = hesson_fit;

if strcmp(fluid,'N2O')
    PDT = load('N2O_PDT_table','PDT');
    PDT = PDT.PDT;
    bubble_rise_velocity_fit = load('N2O_bubble_rise_velocity_fit');
elseif strcmp(fluid,'CO2')
    PDT = load('CO2_PDT_table','PDT');
    PDT = PDT.PDT;
    bubble_rise_velocity_fit = load('CO2_bubble_rise_velocity_fit');
else
    error('fluid string incorrect. try N2O or CO2')
end
constants.bubble_rise_velocity_fit = bubble_rise_velocity_fit.fittedmodel;

% need this code for using qinterp2:
Tvec_table = PDT.T;
Pvec_table = PDT.P;

[Tgrid_table, Pgrid_table] = meshgrid(Tvec_table, Pvec_table);
PDT.T = Tgrid_table;
PDT.P = Pgrid_table;

%% initialize things

[adaptive, a, b, c, bs, s] = butcher_tableau(ode_solver);

constants.error_detected = 0;

k = zeros(N_dim,s);
ti = 0;
n = 1;              % [] counter

t = 0;

if strcmp(ADQMOM, 'off')
    p = 1;
    constants.ADQMOM_p = p;
end

if p == 1
    constants.V_moment_index = 4;
elseif rem(p*3,1) == 0
    constants.V_moment_index = 1 + 3*p;
end

V_moment_index = constants.V_moment_index;

alpha_mom = 4/3 * pi * IC_moments(V_moment_index);

IC_moments = IC_moments * alpha_ic/alpha_mom;

% [r_ic, w_ic] = PD_method(IC_moments);

[r_ic, w_ic] = PD_method_alternative(IC_moments,p);

% r_ic = logspace(-12,-9,N_ab);
%
% w_ic = 1e20*ones(N_ab,1);

if ~isreal(r_ic)
    error('imaginary IC')
end

r_ic = repmat(r_ic(:), N_nodes, 1);
w_ic = repmat(w_ic(:), N_nodes, 1);

DQMOM_IC = [w_ic; r_ic.*w_ic];


% reshape so each ROW is one node
w_ic = reshape(w_ic, N_ab, N_nodes)';

% reshape so each ROW is one node
r_ic = reshape(r_ic, N_ab, N_nodes)';

for i = 1:2*N_ab
    mom(:, i) = sum( r_ic.^((i-1)/p) .* w_ic, 2 );
end

% net volume of bubbles, per unit volume of liquid


V_bubi = 4/3*pi*mom(:,V_moment_index);


% [P] = refpropm('P','T',Ti-0.01,'Q',0.5,fluid);
%
% % [rho_tg, P] = refpropm('-P','T',Ti+0.05,'Q',0.5,fluid);
% % [rho_l] = refpropm('+','T',Ti,'Q',0.5,fluid);
% % [u_tg] = refpropm('U', 'T', Ti +0.05,'Q', 1, fluid);
% P = P*1e3;
%
% rho_l = qinterp2(PDT.T, PDT.P, PDT.D_liq, Ti, P/1e3);
%
% [~, rho_tg, ~, u_tg] = fits_for_getting_P(P, fluid);
%
% rho_tg_sat = rho_tg;

% P = 1e3*refpropm('P','T',Ti,'Q',0,fluid);
%
% [rho_tg_l, rho_tg_v, ~, u_tg_v] = n2o_fits_for_getting_P(P);
% u_tg = u_tg_v;
% % rho_l = rho_tg_l;
% rho_tg = rho_tg_v;
%
% rho_l = qinterp2(PDT.T, PDT.P, PDT.D_liq, Ti, P/1e3);


T_s = Ti;

V_l = fill_level*V_tank;

node_level = get_node_levels(V_l, V_bubi, V_node, V_l/V_node);
guesses.node_level = node_level;

V_bub = sum_over_nodes(V_bubi, node_level, V_node);


% V_bub = V_node/2*V_bubi(1)*node_level(1) + sum(node_level(2:end-1).*V_bubi(2:end-1)*V_node) + V_node/2*V_bubi(end)*node_level(end);

V_l_star = V_bub + V_l;

% V_l_star = fill_level.*V_tank;
% V_bub = V_l_star .* V_bubi;

% V_tg = (1 - fill_level).*(V_tank - V_bub);

V_tg = V_tank - V_l_star;

% V_l = V_l_star.*(1 - V_bubi);

guesses.rho_l = refpropm('D','T',Ti,'Q',0,fluid);

Pi = converge_on_IC(Ti, V_tg, V_l, ...
    V_tank, V_node, V_bubi, PDT, guesses, constants, fluid);

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

P = Pi;
rho_tg_sat = rho_tg;

m_tg = rho_tg * V_tg;
m_l = rho_l * V_l;

U_tg = m_tg * u_tg;
T_l = Ti;

x_tg = 1;
% T_tg = Ti;

T_lw = Ti*ones(N_rw,1);
T_gw = T_lw;

y_i_basic = [m_tg; U_tg; T_gw;  m_l; Ti; T_lw];
% y((N_dim - N_mom + 1):(N_dim), 1) = DQMOM_IC(:);
y_i_DQMOM = DQMOM_IC(:);

y = [y_i_basic; y_i_DQMOM];
% y(7:(6+length(DQMOM_IC)), 1) = DQMOM_IC(:);


% 1 = m_tg
% 2 = U_tg
% 3 = T_gw
% 4 = m_l
% 5 = T_l
% 6 = T_lw
% 7:(N + 6) = weights
% N+7 : 2N+6 = weighted abscissas

[K_b, N_A, h_planck] = universal_constants('boltzmann', 'avagadro', 'planck');

[P_cr, T_cr] = refpropm('PT', 'C', 0, '', 0, fluid);
P_cr = P_cr*1e3;

constants.D = D;
constants.t_w = t_w;
constants.rho_w = rho_w;
constants.cv_w = cv_w;
constants.Cd = Cd;
constants.A_inj = A_inj;
constants.Po = Po;
constants.T_air = T_air;
constants.V_tank = V_tank;
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
V_bub = 0;
gas_holdup = 0;

deltaT_sup = 0;
Pdot = 0;
rhodot_l = 0;
rhodot_tg = 0;

net_nuc = zeros(N_nodes,1);

% I don't seem to be using this term anymore so just leave it
derivatives = zeros(5,1);
guesses.rhodot_tg_sat = 0;
% initialize f
f = zeros(N_dim,1);

constants.min_flag = 0;
constants.peak_flag = 0;

f_cutoff_norm = 1e-2; % supposed to be filter cutoff f / (1/2 sample f)
% so it's basically = 2*dt*f_cutoff
filter_order = 3;
[b_filter,a_filter] = butter(filter_order, f_cutoff_norm, 'low');


%% ODE solver

% begin looping
while running == 1;
    %     figure(1)
    %     hold on
    %     plot(t,P,'bo')
    %     if ~stop.requested
    
    starti = max([n-3, 1]);
    
    Pdot = 0.5*Pdot + 0.5*bdiff(P,starti,n,t,adaptive);
    %     rhodot_l = 0.1*rhodot_l + 0.9*bdiff(rho_l,starti,n,t,adaptive);
    %     rhodot_tg = 0.1*rhodot_tg + 0.9*bdiff(rho_tg,starti,n,t,adaptive);
    Vdot_l(n+1) = V_tank*bdiff(V_l/V_tank,starti,n,t,adaptive);
    %     Vdot_tg(n+1) = V_tank*bdiff(V_tg/V_tank,starti,n,t,adaptive);
    
    if n > 10
        
        rho_tg_sat_filtered = filtfilt(b_filter, a_filter, rho_tg_sat);
        
        guesses.rhodot_tg_sat = bdiff(rho_tg_sat_filtered, starti, n, t, adaptive);
        
    else
        
        guesses.rhodot_tg_sat = 0;
    end
    
    Vdot_bub(n+1) = bdiff(V_bub, starti, n, t, adaptive);
    
    if h > 5*h_min
        
        guesses.Vdot_l = 0.5*Vdot_l(n+1) + 0.5*guesses.Vdot_l;
        
    else
        
        guesses.Vdot_l = 0;
    end
    
    constants.outerloop_superheat = deltaT_sup(n);
    
    % check for min
    if (t(n) > 0.1) && (Pdot > 10)
        %         Pdot
        if constants.min_flag == 0;
            constants.min_flag = 1;
            disp('P min')
            t_min = t(n);
            n_min = n;
            %             running = 0;
        end
    end
    
    % check for peak
    if (constants.min_flag == 1) && (constants.peak_flag == 0)
        if ((Pdot) < 0) && (t(n) > 1.25*t_min)
            constants.peak_flag = 1;
            %             running = 0;
            disp('P peak')
            n_peak = n;
            t_peak = t(n);
        end
    end
    dt = t(n) - t(max([1, n-1]));
    constants.t = t(n);
    constants.dt = dt;
    
    %     if peak_flag == 1
    %         if t(n) > t_peak*1.5
    %             %             running = 0;
    %         end
    %     end
    %
    %         if t(n) > 0.05
    %             running = 0;
    %         end
    
    
    if P(end) < 2e6
        running = 0;
        stop_reason = 'P got below 2 MPa';
    end
    
    
    mdot_l = f(4);
    %     rhodot_l = mdot_l/V_l(n) - m_l(n)/V_l(n)^2*Vdot_l(n);
    
    mdot_tg = f(1);
    %     rhodot_tg = mdot_tg/V_tg(n) - m_tg(n)/V_tg(n)^2*Vdot_tg(n);
    
    
    %     derivatives = zeros(5,1);
    
    
    % 1 = m_tg
    % 2 = U_tg
    % 3 = T_gw
    % 4 = m_l
    % 5 = T_l
    % 6 = T_lw
    % 7:(N + 6) = weights
    % N+7 : 2N+6 = weighted abscissas
    
    % if I'm just playing around, print status at each step
    if nargin == 0
        
        fprintf(['t = %#4.4g, dt = %#4.4g, P = %#4.4g, alpha = %#4.4g, dT_sup = %#6.6g,' ...
            ' V_bub = %#4.4g, T_l = %#4.4g, m_l/m_li = %#4.4g,'...
            ' m_tg/m_tgi -1 = %#4.4g, fill_level%% = %#4.4g, Vdot_l = %#4.4g\n'],...
            t(n), t(n) - t(max([1, n-1])), P(n)/1e6, gas_holdup(n), deltaT_sup(n), ...
            V_bub(n), T_l(n), m_l(n)/m_l(1), m_tg(n)/m_tg(1)-1, 100*fill_level(n), Vdot_l(n+1));
        
        %                 fprintf('y:\t \t \t \t');
        %                 fprintf('%3.3e\t', y(1:6,end))
        %                 fprintf('|')
        %                 fprintf('%3.3e\t', y(7:end,end))
        %
        %
        %         if n > 1
        %         fprintf('\nrelative change in last step:\t')
        %         fprintf('%3.3e\t', (y(1:6,end) - y(1:6,end-1))./(y(1:6,end-1)))
        %         fprintf('|')
        %         fprintf('%3.3e\t', (y(7:end,end) - y(7:end,end-1))./(y(7:end,end-1)))
        %                 fprintf('\n')
        %
        %         end
        
    end
    
    % if we're not on the first step and error is plenty small, increase
    % step size
    if  n > 1 && adaptive == 1
        
        % if error is < min_error
        if max(abs_err/abs_tol,rel_err/rel_tol) < min_error
            
            % make h bigger
            h = min(4*h,h_max);
            
        end
        
        % also check if we're close to LRO
        
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
        % maybe should modify this for going the other way too?
        
        %         if dT_superheat(n) < 0
        
        % slope of dT_superheat curve
        dTs_slope = bdiff(deltaT_sup,starti,n,t,adaptive);
        
        % projected t_sup (t when dT_sup = 0)
        t_sup = -deltaT_sup(n)/dTs_slope + t(n);
        
        
        h_sup = t_sup - t(n); % distance to t_sup
        
        % if the step we're about to take is >3/4 the distance to the
        % crossover point and the distance is bigger than the tolerance
        if (dTs_slope > 0 && deltaT_sup(n) < 0) || (dTs_slope < 0 && deltaT_sup(n) > 0 )
            
            if (h > 0.5*h_sup && h_sup > dT_sup_tol) && (h_sup > 0);
                
                % set h to 1/2 the distance to LRO (ie refine)
                h = 0.5*h_sup;
                disp('refining based on superheat being near 0')
                
            end
            
        end
        
        
    end
    
    error_OK = 0;
    
    while error_OK == 0
        % solving differential equations
        % i = counter for
        
        constants.h = h;
        
        error_flag = 0;
        
        for i = 1:s
            % s = number of stages in the scheme
            
            
            if i == 1
                
                % f = f( t(n) , y(n) )
                
                if n > 1
                    % we calculated it at the end of last time step
                    f = f_np1;
                else
                    % we haven't had a time step yet!
                    constants.step = 1;
                    f = diffeqns(y(:,n), constants, guesses, PDT);
                end
            else
                
                % f for k(2) = f( t(n) + c(2)*h , y(n) + a(2,1)*k(1) )
                % f for k(3) = f( t(n) + c(3)*h , y(n) + a(3,1)*k(1) + a(3,2)*k(2) )
                % and so on
                constants.step = i;
                a_k_term = sum( (ones(N_dim,1)*a(i,1:i-1)).*k(:,1:i-1) ,2 );
                
                y_new = y(:,n) + a_k_term;
                
                variables = unpack_y(y_new, constants);
                
                T_l_new = variables.T_l;
                g_q_new = variables.g_q;
                w_q_new = variables.w_q;
                
                r_q_new = g_q_new./w_q_new;
                
                V_bubi_new = 4/3*pi*sum( r_q_new.^((4-1)/p) .* w_q_new, 2 );
                
                % check to see if we're taking a step that'll cause errors
                error_conditions = [isnan(sum(y_new(:))), (T_l_new < 220), (T_l_new > 305), ...
                    (T_l_new > (T_l(n) + 0.25)),  (T_l_new < (T_l(n) - 5)),...
                    (max(r_q_new(:)) > 0.2), (max(V_bubi_new) > 1)];
                % (min(y_new) < 0),
                
                if sum(error_conditions) > 0
                    error_conditions(:)'
                    disp('we''re taking a bad step')
                    f = ones(size(f));
                    error_flag = 1;
                    break
                else
                    
                    [ind_neg, min_val] = min(y_new);
                    if min_val < 0
                        fprintf('negative part of y. index: %0.d\n', ind_neg)
                        %                     keyboard
                    end
                    
                    f = diffeqns(y_new, ...
                        constants, guesses, PDT);
                    
                end
            end
            
            k(:,i) = f*h;
            
        end
        
        %         k1 = h*diffeqns(y(:,n));
        %         k2 = h*diffeqns(y(:,n) + a(2,1)*k1);
        %         k3 = h*diffeqns(y(:,n) + a(3,1)*k1 + a(3,2)*k2);
        %         k4 = h*diffeqns(y(:,n) + a(4,1)*k1 + a(4,2)*k2 + a(4,3)*k3);
        %         k5 = h*diffeqns(y(:,n) + a(5,1)*k1 + a(5,2)*k2 + a(5,3)*k3 + a(5,4)*k4);
        %         k6 = h*diffeqns(y(:,n) + a(6,1)*k1 + a(6,2)*k2 + a(6,3)*k3 + a(6,4)*k4 + a(6,5)*k5);
        %
        %         k = [k1, k2, k3, k4, k5, k6];
        
        y(:,n+1) = y(:,n) + (k*b);
        
        if error_flag == 0
            % need to check the new y to make sure it's OK too
            
            variables = unpack_y(y(:,n+1), constants);
            
            T_l_new = variables.T_l;
            g_q_new = variables.g_q;
            w_q_new = variables.w_q;
            
            r_q_new = g_q_new./w_q_new;
            
            % check to see if we're taking a step that'll cause errors
            error_conditions = [isnan(sum(y_new(:))), (T_l_new < 220), (T_l_new > 305), ...
                (T_l_new > (T_l(n) + 0.25)),  (T_l_new < (T_l(n) - 5)),...
                (max(r_q_new(:)) > 0.2), (max(V_bubi_new) > 1)];
            % (min(y_new) < 0),
            
            if sum(error_conditions) > 0
                disp('we''re taking a bad step')
                f = ones(size(f));
                error_flag = 1;
            end
        end
        
        if adaptive == 1
            % using adaptive scheme, need to check error
            
            err = k*(b - bs);   % absolute error (diff. between 5th and 4th order estimates of y(n+1) - y(n))
            
            %             rel_err = abs(err./( y(:,n) + 1e-6));  % relative error
            
            %             rel_err = abs(err./( mean([y(:,n) y(:,n+1)],2) + 1e-6));  % relative error
            
            
            for j = 1:N_dim
                if abs(y(j,n)) > 1e-6
                    
                    rel_err(j) = abs(err(j))./( abs( mean( y(j,n:n+1))) + 1e-6);  % relative error
                else
                    rel_err(j) = abs(err(j));
                end
            end
            
            %             rel_err = rel_err(1:6); % remove the bubble distribution terms
            
            [rel_err, ind_max_rel_err(n+1)] = max(rel_err(isfinite(rel_err)));  % fix rel_err to the maximum finite value of rel_err
            
            
            
            abs_err = abs(err);
            
            %             abs_err = abs_err(1:6); % remove the bubble distribution terms
            
            abs_err = max(abs_err(isfinite(abs_err)));  % do the same for abs_err
            
            % check for possible problems: isempty statements are in case
            % abs and rel err are both full of non-finite values
            % isnan checks for nan's
            % isreal checks for imaginary numbers
            error_conditions2 = isempty(rel_err) + ...
                isempty(abs_err) +  ...
                isnan(sum(err)) + ...
                ~isreal(sum(y(:,n+1))) + ...
                error_flag;
            
            % if any of those fail, set rel_err large so that the step gets
            % recomuputed
            if error_conditions2 > 0
                rel_err = 1;
            end
            
            if ( rel_err < rel_tol && abs_err < abs_tol) || (h < 1.25*h_min)
                % meeting the error requirement or step size is too
                % small already
                error_OK = 1;
                
                if ((n > 1) && ((h_LRO < LRO_tol) && (h_LRO > 0))) && (fill_level(n) < 0.01)
                    % distance to LRO is less than LRO_tol
                    running = 0;
                    %                     disp('reached LRO')
                end
                
                if h < 2*h_min
                    fprintf('h got too small. exceeded tolerance by %6.4g%%\n',100*rel_err/rel_tol);
                    running = 0;
                end
                
            else
                
                %                 fprintf(['max rel err = %6.4g, ind of max rel err = %6.4g\n'...
                %                     'err(ind_max_rel_err) = %8.6g, y(ind_max_rel_err,n+1) = '...
                %                     '%8.6g, y(ind_max_rel_err,n) = %8.6g\n'], ...
                %                     rel_err, ind_max_rel_err, err(ind_max_rel_err), ...
                %                     y(ind_max_rel_err, (n+1):-1:n))
                
                % not meeting error requirements
                % sh is used to update h, h = sh*h
                
                if rel_err == 0 || abs_err == 0
                    % something odd happened, so reduce step size a lot
                    
                    sh = 0.1;
                    
                else
                    
                    if rel_err/rel_tol > abs_err/abs_tol
                        % if relative error is a bigger problem than
                        % absolute, update step size based on relative
                        % error. Else, use absolute
                        
                        sh = 0.84*( rel_tol*h / (2*rel_err) )^(1/4);
                    else
                        sh = 0.84*( abs_tol*h / (2*abs_err) )^(1/4);
                    end
                end
                
                if sh < 0.1
                    % if it looks like the step size would be reduced too
                    % much, only reduce it by 1/10
                    sh = 0.1;
                elseif sh > 4.0
                    % similarly if it's too big, only make it 4x bigger
                    sh = 4.0;
                end
                
                % update step size
                h = h*sh;
                
                % minimum step size set by computer's precision
                h_min = 16*eps(t(n));
                
                % self explanatory I think
                if h > h_max
                    h = h_max;
                elseif h < h_min
                    h = h_min;
                end
                
            end
            
        else
            % not using adaptive scheme, don't need to check error
            error_OK = 1;
        end
        
    end
    
    constants.step = 1;
    [f_np1, debug_data] = diffeqns(y(:,n+1), constants, guesses, PDT);
    
    m_tg(n+1) = debug_data.m_tg;
    U_tg(n+1) = debug_data.U_tg;
    m_l(n+1) = debug_data.m_l;
    T_l(n+1) = debug_data.T_l;
    g_q(:,:,n+1) = debug_data.g_q;
    r_q(:,:,n+1) = debug_data.r_q;
    mom(:,:,n+1) = debug_data.mom;
    w_q(:,:,n+1) = debug_data.w_q;
    V_bubi(:,n+1) = debug_data.V_bubi;
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
    u_rise(:,:,n+1) = debug_data.u_rise;
    N_bubi(:,n+1) = debug_data.N_bubi;
    N_full(n+1) = debug_data.N_full;
    
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
    if max(r_q(:)) > 0.99
        disp('r_q is getting too big')
%         keyboard
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
            x_node = [1:length(V_bubi)];
            plot(x_node, N_bubi, 'k-s', x_node, V_bubi, 'k-o', x_node, (6*V_bubi./N_bubi/pi).^(1/3),'k-*')
            set(gca,'yscale','log')
            legend('number density','volume density','SMD')
            
            %             % wall temperature profiles
            %             plot([1:length(T_lw)], T_lw,'k-',[1:length(T_gw)],T_gw,'b')
            
            
            
            
            clock_plot = clock;
        end
    end
    
    
    
    guesses.P = P(n+1);
    guesses.rho_tg = rho_tg(n+1);
    guesses.rho_l = rho_l(n+1);
    guesses.node_level = node_level;
    
    t(n+1) = t(n) + h;
    
    n = n + 1;
    
    if deltaT_sup(n) > 7.5
        running = 0;
        stop_reason = 'superheat got bigger than 7.5';
    end
    
    if (sum(abs(imag(y(:,n)))) > 0) || (sum(isnan(y(:,n))) > 0)
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
    
    
    if save_periodically
        
        clock_now = clock;
        
        time_since_save = etime(clock_now, clock_save);
        if time_since_save >= t_save;
            if save_stuff == 1
                disp('saving')
                save(save_filename,'-v7.3')
            end
            clock_save = clock;
        end
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
%% plotting and output

if save_stuff == 1
    cd(results_save_dir)
    
    if save_parameters_only == 1
        
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
        file_exists = 1;
        file_num = 1;
        while file_exists
            save_filename_numbered = [save_filename num2str(file_num) '.mat'];
            if ~exist(save_filename_numbered,'file')
                save(save_filename_numbered,'-v7.3')
                file_exists = 0;
                cd(current_dir);
            end
            file_num = file_num+1;
        end
        
    end
    
end



if plot_stuff == 1
    
    m_out = cumtrapz(t, mdot_out_liq + mdot_out_vap);
    mh_out = cumtrapz(t, (mdot_out_liq.*h_l + mdot_out_vap.*h_tg_sat));
    
    [P_exp, T_lw_out_exp, LL_exp] = load_experimental_data(t, specified_case);
    
    T_l = y(4+N_rw,:);
    m_l = y(3+N_rw,:);
    
    m_tg = y(1,:);
    
    T_lw_in = y(5+N_rw,:);
    T_gw_in = y(3,:);
    T_lw_out = y(4+2*N_rw,:);
    
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
    plot(t, ind_max_rel_err)
    xlabel('Time [s]')
    ylabel('index []')
    title('index of max relative error')
    
    
    
end

beep

toc

%% differential equations
function varargout = diffeqns(y, constants, guesses, PDT)

% retrieve constants
D = constants.D;
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
V_node = 0.25*pi*D^2*L_node;
N_nodes = constants.N_nodes;
fluid = constants.fluid;
N_ab = constants.N_ab; % number of abscissas
N_mom = 2*N_ab;
hesson_fit = constants.hesson_fit;
bubble_rise_velocity_fit = constants.bubble_rise_velocity_fit;

p = constants.ADQMOM_p;

% retrieve derivatives calculated with backwards differencing
% Pdot = derivatives(1);
% rhodot_l = derivatives(2);
% rhodot_tg = derivatives(3);
% Vdot_l = derivatives(4);
% Vdot_tg = derivatives(5);

% check for negative values
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

variables = unpack_y(y, constants);

% if constants.t >3.738
%     keyboard
% end

m_tg = variables.m_tg;
U_tg = variables.U_tg;
T_gw = variables.T_gw;
m_l = variables.m_l;
T_l = variables.T_l;
T_lw = variables.T_lw;

g_q = variables.g_q;
w_q = variables.w_q;

if isnan(sum(y)) || ~isreal(sum(y))
    disp('problem: nans or imaginary y')
end

% get the abscissas from the weighted abscissas
r_q = g_q./w_q;

if min(r_q) < 0
    disp('negative abscissa')
%     keyboard
end

% get the moments from gauss quadrature
for i = 1:2*N_ab
    mom(:, i) = sum( r_q.^((i-1)/p) .* w_q, 2 );
end

mom(mom<0) = 0;

N_bubi = mom(:,1);

if ~isreal(mom)
    disp('imaginary moments')
end

V_moment_index = constants.V_moment_index;

% bubble volume per unit volume of liquid/bubble mixture (hence the i)
% (can also view this as the vapor volume fraction aka gas holdup)
V_bubi = 4/3*pi*mom(:,V_moment_index);

if sum(imag(V_bubi)) > 0
    disp('V_bubi went imaginary')
end

if max(V_bubi) > 1
    disp('V_bubi > 1!!')
%     keyboard
    V_bubi(V_bubi>1) = 0.99;
end

% get system pressure
% (assumes pressure is same throughout tank, with no gravity head)
P = get_P_from_mU_mT(m_tg, U_tg, m_l, T_l, V_tank, V_node, ...
    V_bubi, fluid, PDT, constants, guesses);

if (P == pi) || isnan(P)
    disp('P error')
    constants.error_detected = 1;
    P = guesses.P;
%     keyboard
end

% get density of liquid and ullage based on temperature and pressure


if strcmp(constants.property_source,'PDT')
    
    rho_l = qinterp2(PDT.T, PDT.P, PDT.D_liq, T_l, P/1e3);
    
    % if rho_l is NaN, it means we went outside the bounds of PDT, so
    % instead extrapolate it using interp2 (slower than qinterp2)
    if isnan(rho_l)
        rho_l = interp2(PDT.T, PDT.P, PDT.D_liq, T_l, P/1e3, 'spline');
    end
    
    [rho_tg_l, rho_tg_v, u_tg_l, u_tg_v] = fits_for_getting_P(P, fluid);
    
    
elseif strcmp(constants.property_source,'refprop')
    
    rho_l = get_D_from_TP(T_l, P, guesses, constants, fluid);
    
    [rho_tg_l, rho_tg_v, u_tg_v] = refpropm('+-U','P',P/1e3,'Q',1,fluid);
    u_tg_l = refpropm('U','P',P/1e3,'Q',0,fluid);
    
    
end


% rho_tg = qqinterp2(PDT.T, PDT.P, PDT.D_vap, T_tg, P, 'linear');

% [rho_tg, T_tg] = refpropm('DT', 'P', P/1e3, 'U', U_tg/m_tg, fluid);

% get saturation properties for ullage
% at some point should include a switch here to take into account times
% when the ullage is just superheated vapor (not saturated, not metastable)
% [rho_tg_l, rho_tg_v, u_tg_l, u_tg_v] = fits_for_getting_P(P, fluid);

u_tg = U_tg/m_tg;
x_tg = (u_tg - u_tg_l)/(u_tg_v - u_tg_l);
alpha = 1/( 1 + rho_tg_v/rho_tg_l * (1 - x_tg)/x_tg );
rho_tg = alpha*rho_tg_v + (1 - alpha)*rho_tg_l;

% if constants.step == 1
%     fprintf('1 - x_tg = %0.3g\n',1 - x_tg)
% end

% if x_tg > 1
% %     keyboard
%     disp('x_tg > 1')
%     x_tg = 1;
% end

% total liquid volume (not including bubbles)
V_l = m_l./rho_l;

% V_l_star = volume of liquid and bubbles
% = V_l + V_bub
% V_bub = V_bubi*(V_l + V_bub)
% V_l_star = V_l./(1 - V_bubi);

% I've kind of changed the meaning so that V_l_star is in fact the TOTAL
% volume of the liquid+bubbles

% how full each node is. will be 1 for all nodes below liquid level, 0 for
% all nodes above, and [0,1] for the node in which the liquid level
% currently sits
node_level = get_node_levels(V_l, V_bubi, V_node, guesses.node_level);
% V_bub = sum(node_level.*V_bubi*V_node);
V_bub = sum_over_nodes(V_bubi, node_level, V_node);

V_l_star = V_l + V_bub;

% number of nodes still in the liquid (completely)
N_full = sum(node_level == 1);

V_tg = m_tg/rho_tg;

% liquid properties

[h_l, dh_drho_l, drho_dP_l, u_l, Cv_l, dP_dT_l, ...
    k_l, Cp_l, s_l, MW, mu_l] = ...
    refpropm('H!RUO#LCSMV','T',T_l,'D&',rho_l,fluid);



[P_sat, s_liq_sat, h_liq_sat] = refpropm('PSH', 'T', T_l, 'Q', 0, fluid);

dP_drho_l = 1e3./drho_dP_l;
dP_dT_l = dP_dT_l*1e3;
alpha_l = k_l./(rho_l .* Cp_l);
P_sat = 1e3*P_sat;

nu_l = mu_l/rho_l;

Pr_l = nu_l/alpha_l;

% a lot of these properties aren't needed until later, but by calculating
% them here I can reduce the number of refprop calls. Note that the next
% two refprop calls are for ullage properties (liquid and vapor parts)

% properties needed for Vdot calculation (liquid)
[u_tg_l_sat, rho_tg_l, dP_dT_tg_sat, drho_dP_T, drho_dT_P, dh_dT_P, dh_dP_T...
    ,sigma, h_l_sat, T_tg] = refpropm('UDERW(*IHT', 'P', P/1e3, 'Q', 0, fluid);
dP_dT_tg_sat = dP_dT_tg_sat * 1e3;
drho_dP_T = drho_dP_T * 1e-3;
dh_dP_T = dh_dP_T * 1e-3;

du_dT_P = dh_dT_P + P/rho_tg_l^2 * drho_dT_P;
du_dP_T = dh_dP_T + 1/rho_tg_l + P/rho_tg_l^2 * drho_dP_T;

du_dT_sat_tg_l = du_dT_P + du_dP_T * dP_dT_tg_sat;
drho_dT_l_sat = drho_dT_P + drho_dP_T * dP_dT_tg_sat;
% drho_dP_sat = drho_dP_T + drho_dT_P / dP_dT_sat;


% need to calculate a bunch of partial derivatives for the Vdot function
% extras are needed for the saturated ullage because the 2 phases adds more
% terms.

% properties needed for Vdot calculation (vapor)
[u_tg_v_sat, rho_tg_v, drho_dP_T, drho_dT_P, dh_dT_P, dh_dP_T,...
    T_s, h_tg_sat, s_tg_sat, Cp_tg, mu_tg] = ...
    refpropm('UDRW(*THSCV', 'P', P/1e3, 'Q', 1, fluid);
drho_dP_T = drho_dP_T * 1e-3;
dh_dP_T = dh_dP_T * 1e-3;

du_dT_P = dh_dT_P + P/rho_tg_v^2 * drho_dT_P;
du_dP_T = dh_dP_T + 1/rho_tg_v + P/rho_tg_v^2 * drho_dP_T;

du_dT_sat_tg_v = du_dT_P + du_dP_T * dP_dT_tg_sat;
drho_dT_v_sat = drho_dT_P + drho_dP_T * dP_dT_tg_sat;
% drho_dP_sat_tg = drho_dP_T + drho_dT_P / dP_dT_sat;

drho_dx_P_tg = -rho_tg^2 *(1/rho_tg_v - 1/rho_tg_l);
drho_dP_x_tg = (1/dP_dT_tg_sat) * rho_tg^2 * ( x_tg/rho_tg_v^2 * drho_dT_v_sat + ...
    (1-x_tg)/rho_tg_l^2 * drho_dT_l_sat );

rho_tg_sat = rho_tg_v;
% % temp of saturated surface based on pressure (and h of sat. vapor)
% [T_s, h_tg_sat, rho_tg_sat] = refpropm('THD','P',P/1e3,'Q',1,fluid);
%
% % saturated liquid enthalpy at P
% [sigma, h_l_sat] = refpropm('IH','P',P/1e3,'Q',0,fluid);

% heat of vaporization (at saturation based on P)
h_lv = h_tg_sat - h_l_sat;

% bubble calculations

% superheat = T - T_sat
deltaT_sup = T_l - T_s;

% if abs(constants.outerloop_superheat - deltaT_sup)/deltaT_sup > 5e-2
%     fprintf('superheat of outer loop and diff eqns differs. \n outer: %4.4g, inner: %4.4g\n', constants.outerloop_superheat, deltaT_sup)
% end

% mass flow rate out via injector

% passing 2-phase properties
% x = rho_tg_sat/rho_l / (1/V_bubi(1) + rho_tg_sat/rho_l - 1);
x_out = constants.C_x_inj * V_bubi(1)/(V_bubi(1) + rho_l/rho_tg_sat*(1 - V_bubi(1)));

% if x_inj < 0.05
%     x_inj = x_inj + 0.05;
% end

rho_liq_mix = rho_l*(1-V_bubi(1)) + rho_tg_sat*V_bubi(1);
s_liq_mix = s_l * (1-x_out) + s_tg_sat * x_out;
h_liq_mix = h_l * (1-x_out) + h_tg_sat * x_out;

% mdot_out_mix = A_inj*Cd*injector_flow(x_inj, P, hesson_fit);

mdot_out_mix = A_inj*Cd*injector_flow(x_out, P, hesson_fit, rho_liq_mix, h_liq_mix, fluid, constants);
mdot_out_liq = (1 - x_out)*mdot_out_mix;
mdot_out_vap = x_out*mdot_out_mix;
% passing straight liquid properties
% from bottom node
% mdot_out = A_inj*Cd*injector_flow(Po, P, T_l, rho_l, P_sat, s_l, h_l);

% bulk flow velocity out the bottom
u_bulk = mdot_out_mix / rho_l / (0.25 * pi * D^2);
u_bulk = 0;

Mo = g*mu_l^4*(rho_l - rho_tg_sat)/(rho_l^2*sigma^3);

Ja = Cp_l * rho_l * deltaT_sup/(rho_tg_v * h_lv);

% rise velocity due to buoyancy
% u_rise = sqrt( 2.14*sigma./( rho_l*(2*abs(r_q) ) + 0.505*g*( 2*abs(r_q) ) ) );
% u_rise = 4 * 0.71*sqrt(g*2*abs(r_q));
% u_rise = 2*g*rho_l*r_q.^2./(9*sigma);


% using fits compiled from complicated correlations
% u_rise = constants.C_u_rise * bubble_rise_velocity(r_q, T_tg);
% u_rise = constants.C_u_rise * exp( bubble_rise_velocity_fit(log(r_q), T_tg*ones(size(r_q)) ));
%
% % have to do this because the velocity fit is bad in this region
% u_rise(r_q < 1e-6) = 0;

u_rise = zeros(size(r_q));

delta_rho = rho_l - rho_tg_sat;

for i = 1:N_full+1
    
    if V_bubi(i) > 0.3
        Eprime = (1 - V_bubi(i)).^2;
        Cd = 8*Eprime/3;
        u_rise(i,:) = sqrt(8/3 * delta_rho * g * r_q(i,:)./(rho_l * Cd));
    else
        
        
        % fit from fan and tsuchiya
        %
        n_FT = 1.6;
        c_FT = 1.2;
        Kbo_FT = 10.2;
        Kb_FT = max([ Kbo_FT*Mo^-0.038, 12]);
        
        De = 2*r_q(i,:);
        
        De_ND = De*sqrt(rho_l*g/sigma);
        
        
        u_rise_ND = ( (Mo^-0.25/Kb_FT * ( delta_rho/rho_l )^(5/4) * De_ND.^2 ).^-n_FT ...
            + ( 2*c_FT./De_ND + ( delta_rho/rho_l ) * De_ND/2 ).^(-n_FT/2) ).^(-1/n_FT);
        
        u_rise_FT = u_rise_ND/(rho_l/(sigma*g))^(1/4);
        
        % u_rise = constants.C_u_rise * u_rise;%.*(( (1 - V_bubi.^(5/3))./(1 - V_bubi).^2) * ones(size(u_rise(1,:))) ).^(-0.5);
        
        alpha_IZ = V_bubi(i)*ones(1,N_ab);
        
        mu_mix = mu_l*(1 - alpha_IZ).^(-2.5*(mu_tg + 0.4*mu_l)/(mu_tg + mu_l));
        
        Re_inf = rho_l * u_rise_FT .* 2 .* r_q(i,:)./mu_l;
        
        u_rise(i,:) = u_rise_FT.* (1 - alpha_IZ).*mu_l./mu_mix .* ...
            (1 + 0.1*Re_inf.^0.75)./(1 + 0.1*Re_inf.^0.75.*( sqrt(1 - alpha_IZ).*mu_l./mu_mix).^(6/7));
        
        Re = rho_l*u_rise(i,:).*2.*r_q(i,:)./mu_mix;
        
        for j = 1:N_ab
            
            if Re(j) > 1e3
                E = ( (1 + 17.67*(1 - V_bubi(i)).^(6/7) )./(18.67*(1 - V_bubi(i)))).^2;
                Eo = g*delta_rho*4*pi*r_q(i,j).^2/sigma;
                Cd = 2/3*E.*Eo;
                u_rise(i,j) = sqrt(8/3 * delta_rho * g * r_q(i,j)./(rho_l * Cd));
                
            end
            
        end
        
    end
    
end


u_rise(N_full+2:end,:) = 0;

if constants.step == 1
    
    u_max = max(u_rise(:));
    CFL = constants.h*u_max/L_node;
    if CFL > 0.5
        fprintf('CFL>0.5, = %0.3g\n', CFL);
    end
end

Re_max = max(max(r_q(1:N_full+1,:).*u_rise(1:N_full+1,:)))*rho_l/mu_l;
Re_min = min(min(r_q(1:N_full+1,:).*u_rise(1:N_full+1,:)))*rho_l/mu_l;

Eo = g*(rho_l - rho_tg_sat)*r_q.^2*4/sigma;

Eo_max = max(max(Eo(1:N_full+1,:)));
Eo_min = min(min(Eo(1:N_full+1,:)));

% if constants.step == 1 && constants.t > 0.001
%     fprintf('Re: %0.2g - %0.2g, Eo: %0.2g - %0.2g, Mo: %0.2g, Ja: %0.2g\n', Re_min, Re_max, Eo_min, Eo_max, Mo, Ja)
% end

% u_rise(r_q < 1e-6) = 0;


if sum(sum(isnan(u_rise(1:N_full+1,:)))) > 0
    disp('nans in u_rise')
%     keyboard
end

duw_dx = zeros(size(w_q));
dug_dx = duw_dx;

u_vec = u_rise - u_bulk;
uw_vec = (u_rise - u_bulk).*w_q;
ug_vec = (u_rise - u_bulk).*g_q;


u_physical = max(abs(u_rise(:) - u_bulk));
CFL = u_physical * constants.h / L_node;

if CFL > 0.9
    disp('CFL is greater than 0.9 (inside diff eqns)')
end

for i = 1:N_full + 1
    % fluxes in and out of node
    
    %         % using finite_diff (various possibilities) - play around w/
    %         forwards, backwards, central, and what order
    %         max_order = 3;
    %         if i > 1
    %
    %
    %             if i == N_full + 1
    %                 for j = 1:N_ab
    %                     duw_dx(i,j) = finite_diff(uw_vec(:,j), 1, N_full+1, i, L_node, 'backwards', max_order);
    %                     dug_dx(i,j) = finite_diff(ug_vec(:,j), 1, N_full+1, i, L_node, 'backwards', max_order);
    %                 end
    %             else
    %
    %                 for j = 1:N_ab
    %                     duw_dx(i,j) = finite_diff(uw_vec(:,j), 1, N_full+1, i, L_node, 'backwards', max_order);
    %                     dug_dx(i,j) = finite_diff(ug_vec(:,j), 1, N_full+1, i, L_node, 'backwards', max_order);
    %                 end
    %             end
    %
    %         else
    %             % bottom node
    %             duw_dx(i,:) = uw_vec(i,:)/L_node;
    %             dug_dx(i,:) = ug_vec(i,:)/L_node;
    %         end
    
    %
    
    % %     1st order upwind
    %     if i > 1 && i < N_full
    % %         interior grid points
    %
    %         for j = 1:N_ab
    %
    %             if u_rise(i,j) - u_bulk > 0
    %                 rising: flux in is from -x
    %                 flux_bot = w_q(i-1,j)*0.5*(u_vec(i,j) + u_vec(i-1,j));
    %                 flux_top = w_q(i,j)*0.5*(u_vec(i+1,j) + u_vec(i,j));
    %                 duw_dx(i,j) = (flux_top - flux_bot)/L_node;
    %
    %                 flux_bot = g_q(i-1,j)*0.5*(u_vec(i,j) + u_vec(i-1,j));
    %                 flux_top = g_q(i,j)*0.5*(u_vec(i+1,j) + u_vec(i,j));
    %                 dug_dx(i,j) = (flux_top - flux_bot)/L_node;
    %             else
    %                 falling: flux in is from +x
    %                 flux_top = w_q(i+1,j)*0.5*(u_vec(i+1,j) + u_vec(i,j));
    %                 flux_bot = w_q(i,j)*0.5*(u_vec(i,j) + u_vec(i-1,j));
    %                 duw_dx(i,j) = (flux_top - flux_bot)/L_node;
    %
    %                 flux_top = g_q(i+1,j)*0.5*(u_vec(i+1,j) + u_vec(i,j));
    %                 flux_bot = g_q(i,j)*0.5*(u_vec(i,j) + u_vec(i-1,j));
    %                 dug_dx(i,j) = (flux_top - flux_bot)/L_node;
    %             end
    %
    %         end
    %
    %     else
    %
    %         if i == 1
    %             bottom point
    %
    %             duw_dx(i,:) = uw_vec(i,:)/L_node;
    %             dug_dx(i,:) = ug_vec(i,:)/L_node;
    %                         duw_dx(i,:) = (uw_vec(i+1,:) - uw_vec(i,:))/L_node;
    %                         dug_dx(i,:) = (ug_vec(i+1,:) - ug_vec(i,:))/L_node;
    %         else
    %             top 2 points
    %                         duw_dx(i,:) = ( uw_vec(i,:) - uw_vec(i-1,:) )/(L_node);
    %                         dug_dx(i,:) = ( ug_vec(i,:) - ug_vec(i-1,:) )/(L_node);
    %             duw_dx(i,:) = zeros(size(uw_vec(i,:)));
    %             dug_dx(i,:) = zeros(size(ug_vec(i,:)));
    %         end
    %     end
    %
    % Roe
    if i > 1 && i < N_full
        % interior grid points
        
        
        for j = 1:N_ab
            % w flux through top
            u_RL = ( sqrt(w_q(i+1,j))*u_vec(i+1,j) + sqrt(w_q(i,j))*u_vec(i,j))...
                /(sqrt(w_q(i+1,j)) + sqrt(w_q(i,j)));
            
            du = u_vec(i+1,j) - u_vec(i,j);
            dw = w_q(i+1,j) - w_q(i,j);
            
            A_ip12 = u_RL;
            
            if abs(u_RL) < du
                phi_RL = (u_RL^2+ du^2)/(2*du);
            else
                phi_RL = abs(u_RL);
            end
            
%             flux_top = 0.5*A_ip12*(w_q(i,j) + w_q(i+1,j)) - 0.5*(phi_RL*dw);
            
            flux_top = 0.5*(uw_vec(i,j) + uw_vec(i+1,j)) - 0.5*(phi_RL*dw);

            
            % w flux through bottom
            u_RL = ( sqrt(w_q(i,j))*u_vec(i,j) + sqrt(w_q(i-1,j))*u_vec(i-1,j))...
                /(sqrt(w_q(i,j)) + sqrt(w_q(i-1,j)));
            
            du = u_vec(i,j) - u_vec(i-1,j);
            dw = w_q(i,j) - w_q(i-1,j);
            
            A_im12 = u_RL;
            
            if abs(u_RL) < du
                phi_RL = (u_RL^2+ du^2)/(2*du);
            else
                phi_RL = abs(u_RL);
            end
            
%             flux_bot = 0.5*A_im12*(w_q(i-1,j) + w_q(i,j)) - 0.5*(phi_RL*dw);
            
            flux_bot = 0.5*(uw_vec(i-1,j) + uw_vec(i,j)) - 0.5*(phi_RL*dw);

            duw_dx(i,j) = (flux_top - flux_bot)/L_node;
            
            
            
            % g flux through top
            u_RL = ( sqrt(g_q(i+1,j))*u_vec(i+1,j) + sqrt(g_q(i,j))*u_vec(i,j))...
                /(sqrt(g_q(i+1,j)) + sqrt(g_q(i,j)));
            
            du = u_vec(i+1,j) - u_vec(i,j);
            dg = g_q(i+1,j) - g_q(i,j);
            
            A_ip12 = u_RL;
            
            if abs(u_RL) < du
                phi_RL = (u_RL^2+ du^2)/(2*du);
            else
                phi_RL = abs(u_RL);
            end
            
%             flux_top = 0.5*A_ip12*(g_q(i,j) + g_q(i+1,j)) - 0.5*(phi_RL*dg);
            
            flux_top = 0.5*(ug_vec(i,j) + ug_vec(i+1,j)) - 0.5*(phi_RL*dg);

            
            % g flux through bottom
            u_RL = ( sqrt(g_q(i,j))*u_vec(i,j) + sqrt(g_q(i-1,j))*u_vec(i-1,j))...
                /(sqrt(g_q(i,j)) + sqrt(g_q(i-1,j)));
            
            du = u_vec(i,j) - u_vec(i-1,j);
            dg = g_q(i,j) - g_q(i-1,j);
            
            A_im12 = u_RL;
            
            if abs(u_RL) < du
                phi_RL = (u_RL^2+ du^2)/(2*du);
            else
                phi_RL = abs(u_RL);
            end
            
%             flux_bot = 0.5*A_im12*(g_q(i-1,j) + g_q(i,j)) - 0.5*(phi_RL*dg);
            
            flux_bot = 0.5*(ug_vec(i-1,j) + ug_vec(i,j)) - 0.5*(phi_RL*dg);

            
            dug_dx(i,j) = (flux_top - flux_bot)/L_node;
        end
        
    else
        
        if i == 1
            % bottom point
            % there's flux out the top and nothing out the bottom
            duw_dx(i,:) = uw_vec(i,:)/L_node;
            dug_dx(i,:) = ug_vec(i,:)/L_node;
            
            %             duw_dx(i,:) = (uw_vec(i+1,:) - uw_vec(i,:))/L_node;
            %             dug_dx(i,:) = (ug_vec(i+1,:) - ug_vec(i,:))/L_node;
        else
            % top 2 points
            %             duw_dx(i,:) = ( uw_vec(i,:) - uw_vec(i-1,:) )/(L_node);
            %             dug_dx(i,:) = ( ug_vec(i,:) - ug_vec(i-1,:) )/(L_node);
            duw_dx(i,:) = zeros(size(uw_vec(i,:)));
            dug_dx(i,:) = zeros(size(ug_vec(i,:)));
        end
    end
    
    
end

if deltaT_sup > 15
    disp('real superheat!')
%     keyboard
end


for i = 1:N_full + 1
    
    % if superheated, then calculate bubble stuff
    
    
    if sum(abs(imag([r_q(i,:); w_q(i,:)]))) > 0
        fprintf('imaginary abscissas or weights. moments:')
        fprintf('%0.6g\t',mom)
        fprintf('\n')
        
    end
    
    %         T_sat = 19.6426*(P/1e3)^0.2499 + 122.3663;
    
    % depth from surface
    depth = (i < N_full)*(N_full - i)*L_node + node_level(N_full+1)*L_node;
    dP_depth = rho_l*g*depth;
    dT_sat_depth = 19.6426*( ((P + dP_depth)/1e3)^0.2499 - (P/1e3)^0.2499 );
    
    deltaT_sup_node = deltaT_sup - dT_sat_depth;
    
    
    if deltaT_sup_node > 1e-4
        % jakob number
        Ja_T = Cp_l * rho_l * C_dTs * deltaT_sup_node/(rho_tg_v * h_lv);
        
        % radius of new bubbles
        r_nuc = C_r_nuc * 2*sigma*T_s/(rho_tg_v * h_lv * C_dTs * deltaT_sup_node);
        
        % bubble radius rate of change
        
        % growth rate for a bubble at rest in an infinite fluid
        % simplified model: plesset & zwick
        % more complicated: scriven
        rdot_rest_plesset = C_rdot * Ja_T^2 * alpha_l ./ r_q(i,:); % bubble at rest
        
        % equation 47 from scriven
        beta_47 = sqrt(0.5*Ja_T./(1 + (Cp_l - Cp_tg)/Cp_l * rho_tg_v/rho_l * Ja_T));
        
        % eq 47 + eq 49 from scriven
        beta_rdot = 0.85*(beta_47 + sqrt(12/pi)*beta_47.^2);
        
        % use r_nuc as a crude approx of r when the bubble first nucleated
        rdot_rest_scriven = 2*beta_rdot.^2*alpha_l./(r_q(i,:));
        
        rdot_rest = rdot_rest_scriven;
        
        % bubble rising in the fluid. from legendre 1998
        rdot_rise = Ja_T * sqrt( 2 * alpha_l * (u_rise(i,:) + 1e-6)./(pi * r_q(i,:) ) ); % rising in the liquid
        
        %                 rdot = max(rdot_rest, rdot_rise);
        rdot = rdot_rest + rdot_rise;
        
        
        
        % length of liquid node volume [m]
        %         L_l = V_l_star(i) / (pi * 0.25 * D^2);
        L_l = node_level(i)*L_node;
        if i == 1 || i == N_nodes
            L_l = L_l/2;
        end
        
        % surface area of node[m^2]
        A_l = pi * D * L_l;% + pi * 0.25 * D^2;
        
        if i == 1
            A_l = A_l + 0.25*pi*D^2;
        end
        
        switch constants.nuc_model
            
            case 'SJ'
                
                % rate of nucleation (mostly from shin & jones, 1993)
                
                % departure diameter [m]
                r_dep1 = 0.5 * 2.97e4 * (P/P_cr)^-1.09 * ( K_b * T_cr / (P_cr * MW) )^(1/3);
                
                Ja_wall = (T_lw(1) - (T_s + dT_sat_depth)) * rho_l * Cp_l/(rho_tg_v*h_lv);
                
                K_1 = (Ja_wall/Pr_l)*( g*rho_l*(rho_l-rho_tg_v)/mu_l^2 * (sigma/(g*(rho_l-rho_tg_v)))^(3/2) )^(-1);
                
                Eo_dep = ( 0.19*( 1.8 + 1e5*K_1)^(2/3) )^2;
                
                r_dep = 0.5 * (Eo_dep*sigma /(g*(rho_l - rho_tg_v)) )^0.5;
                
                %                 if (i == 1 && constants.step == 1) && constants.t > 0.1
                %                     fprintf('r_dep/r_dep1 = %0.4g\n', r_dep/r_dep1)
                %                 end
                
                % correlation from Jensen & Memmel, 1986
                % gives values on the order of 10-20 microns
                
                % non-dimensional cavity size (ie bubble nucleation size)
                r_c_star = r_nuc / r_dep;
                
                % non-dimensional nucleation rate
                N_ns_star = 1e-7 * r_c_star;
                
                % nucleation site density [1/m^2]
                nuc_density = N_ns_star * ( 0.5 / r_dep )^2; % original expression
                
                % other ones I've tried
                %                 nuc_density = N_ns_star * ( 0.5 / r_nuc )^2;
                %                 nuc_density = 1e-4 * ( 0.5 / r_nuc )^2;
                
                
                
                % nucleation frequency [Hz]
                nuc_freq = 1e4 * C_dTs * deltaT_sup_node^n_nuc_freq;
                
                % nucleation rate [Hz]
                nuc_rate = C_nuc_rate * nuc_density * nuc_freq * A_l;
                
            case 'AL'
                
                % alamgir and lienhard, 1981:
                
                B = K_b * T_l / h_planck;
                
                %             phi = 1e-2;
                v_g = 1/rho_tg;
                v_f = 1/rho_l;
                
                J = ( N_A/( 1/rho_l * MW) )^(2/3) * B * ...
                    exp( -16*pi*sigma^3*phi / ( 3*K_b*T_s*(1 - v_f/v_g)^2 * (P_sat - P)^2 ) );
                
                nuc_rate = J * A_l;
                
        end
        
        % homogeneous nucleation rate
        
        rstar = 2*sigma/(P_sat - P);
        
        N_l = rho_l / MW * 1e3*N_A;
        
        m_M = 1e3 * MW / N_A;
        
        J_hom = N_l * rho_l/rho_tg * sqrt(2*sigma/(pi * m_M)) ...
            * exp( - 4*pi*sigma*rstar^2/(3*K_b*T_l));
        
        nuc_rate = nuc_rate + J_hom;
        
        %         if constants.step == 1 && i == 1
        %             fprintf('J_hom/NR = %0.4g\n', J_hom/nuc_rate);
        %         end
        % %
        
        
    else
        if deltaT_sup < -0.5
            disp('real subcooled')
            %             keyboard
        end
        % no superheat -> no bubble growth or nucleation
        r_nuc = 0;
        rdot = 0;
        nuc_rate = 0;
        C_death_rate = 0;
        
    end
    
    % rdot term based on rho_dot
    rdot_rhodot = 0;%-guesses.rhodot_tg_sat/rho_tg_sat * r_q(i,:)/3;
    
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
    
    % nucleation rate per volume
    
    % think this maybe shouldn't have the node level term in it, and
    % maybe not the V node either -> maybe V_l_star?
    % changed my mind -> nuc_rate is based on V_node
    
    spec_nuc_rate = nuc_rate / ( (node_level(i) + 1e-3) * V_node );
    
    if i == 1 || i == N_nodes
        spec_nuc_rate = spec_nuc_rate/(0.5);
    end
    
    %     if (i == 1 && constants.step == 1)
    %         fprintf('spec nuc rate = %0.4g\n growth_int = %0.4g\n',spec_nuc_rate,growth_int_s(V_moment_index)*r_m^3)
    %     end
    
    
    % preallocate
    birth_int_s = zeros(1,N_ab*2);
    
    if deltaT_sup_node > 1e-4
    
    for k = 1:N_ab*2
        % if nucleation happens only at r_nuc
        birth_int_s_delta = (r_nuc/r_m).^((k-1)/p) * spec_nuc_rate;
        
%                 % exponential distribution
%                 r_a = 10*r_nuc;
%         
%                 birth_int_s_exp = (r_a/r_m)^((k-1)/p) * spec_nuc_rate * exp( r_nuc/r_a ) ...
%                     * gamma(1+((k-1)/p)) * gammainc(r_nuc/r_a, 1+((k-1)/p), 'upper');
%         %
%         uniform distribution
%                 dr = 100*r_nuc;
%                 birth_int_s_uni = (1/r_m)^((k-1)/p) * spec_nuc_rate/dr * 1/( (k-1)/p + 1)*...
%                     ( (r_nuc + dr)^( (k-1)/p + 1) - r_nuc^( (k-1)/p + 1) );
        
        
        birth_int_s(k) = birth_int_s_delta;
        
        % now if it's spread out a bit (gaussian)
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
    if strcmp(constants.coalescence_switch,'on')
        
        % most of this is from prince and blanch (1990)
        % buoyancy and turbulence driven
        
        nu_t = 0.0536*D^1.77/rho_l;
        U_max = ( (1 - 0.75*V_bubi(i))/(1 - V_bubi(i)) )*V_bubi(i) * D^2/(48*nu_t);
        mean_shear = 0.53*U_max/(0.5*D);
        
        P2 = P;
        
        liquid_height = V_l_star/(pi*D^2/4);
        
        P1 = P2 + rho_l*g*liquid_height;
        
        L = V_tank/(pi*D^2/4);
        
        Q = V_tank/8; % volumetric gas flow rate. just assumed a constant value here
        
        turb_diss = Q*g * P2*log(P1/P2) / ( pi * (0.5*D)^2 * (P1 - P2) );
        
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
                
                % theta for laminar shear
                qLS = C_coalescence(2) * 4/3*(rbi + rbj).^3*mean_shear;
                
                % rise velocity for bubble i and j
                u_ri = sqrt( (2.14*sigma/(rho_l*dbi)) + 0.505*g*dbi);
                u_rj = sqrt( (2.14*sigma./(rho_l*dbj)) + 0.505*g*dbj);
                
                % collision area
                Sij = pi/4*(rbi + rbj).^2;
                
                % theta for buoyancy
                qB = Sij.*abs(u_ri - u_rj);
                
                % theta for turbulence
                qT = C_coalescence(3) * 0.089*pi*(dbi + dbj).^2 .* turb_diss.^(1/3) .* sqrt(dbi.^(2/3) + dbj.^(2/3));
                
                % bubble radius
                rb_eq = ( (1./rbi + 1./rbj)/2 ).^-1;
                
                % contanct time
                % prince + blanch, 1990, in turn from Levich, 1962
                % note I added the 0.1 factor in front. P+B say this is
                % really just an order of magnitude (if that) estimate
                %                  t_cont = 0.1*rb_eq.^(2/3) ./ turb_diss.^(1/3);
                
                % kamp & chesters, 2001
                rho_c = rho_l;
                C_vm = 0.8;
                t_cont = sqrt( rho_c*C_vm/(3*sigma) * ( 2*dbi.*dbj./(dbi + dbj)).^3 );
                
                % film initial and final thicknesses
                film_i = 1e-4;
                film_f = (C_hamaker * rb_eq./(8*pi*sigma)).^(1/3);
                
                % time required for coalescence
                t_coal = sqrt( rb_eq.^3 * rho_l/(16 * sigma) ) .* log( film_i ./ film_f);
                
                % coalescence kernel
                beta = C_coalescence(1)*(qT + qB + qLS).*exp( - t_coal ./ t_cont);
                
                
                coal_birth_s(k) = coal_birth_s(k) + sum( 0.5 * w_q(i,l) * w_q(i,:) .* ( abs(r_s(l))^3 + abs(r_s).^3 ).^((k-1)/(3*p)) .* beta );
                coal_death_s(k) = coal_death_s(k) + sum( r_s(l)^((k-1)/p) * w_q(i,l) * w_q(i,:) .* beta );
                
            end
        end
        
        
    end
    
    if abs(coal_birth_s(V_moment_index) - coal_death_s(V_moment_index))/coal_birth_s(V_moment_index) > 1e-6
        
        %         if mom(i,V_moment_index) > 1e-6
        disp('coalescence isn''t conserving mass')
%         keyboard
        %         end
    end
    
    dmom_dt_s = birth_int_s(:) + growth_int_s(:) + coal_birth_s(:) - coal_death_s(:);
    
    r_sp = r_s;
    
    for j = 0:(2*N_ab - 1)
        if j == 0
            A1(1,:) = ones(1,N_ab);
            A2(1,:) = zeros(1,N_ab);
        elseif j == 1
            A1(2,:) = (p-1)/p * r_sp.^(1/p);
            A2(2,:) = (1/p) * r_sp.^(1/p - 1);
        else
            A1(j+1,:) = (1 - j/p) * r_sp.^(j/p);
            A2(j+1,:) = (j/p) * r_sp.^(j/p - 1);
        end
    end
    
    A = [A1 A2];
    
    beta_q = dmom_dt_s;
    
    %     alpha_q = A\beta_q;
    alpha_q = linear_equation_solver(A,beta_q);
    
    a_q = alpha_q(1:N_ab);
    b_q_s = alpha_q(N_ab+1:end);
    %     b_q = b_q_s.*r_q(i,:)'; % go from scaled b (ie b*) back to the real thing
    b_q = b_q_s * r_m;
    
    % have to put this in the right place
    % for i = 1, it should be 1, 2, ... N_ab
    % for i = 2, N_ab+1, N+2, ... N_ab+N_ab
    % for i = i, it's (i-1)*N_ab, ... i*N_ab;%
    ind_node = (i-1)*N_ab + [1:N_ab];
    
    
    % only include birth/death/growth (ie 0D)
    %     dw_dt(ind_node) = a_q;
    %     dg_dt(ind_node) = b_q;
    
    % include flux of bubbles in physical space
    dw_dt(ind_node) = a_q - duw_dx(i,:)';
    dg_dt(ind_node) = b_q - dug_dx(i,:)';
    
    % these have units of (m^3/(m^3*s)) ie (volume/time)/(volume of mix)
    birth_term(i) = birth_int_s(V_moment_index)*r_m^(3);
    growth_term = growth_int_s(V_moment_index)*r_m^(3);
    %     death_term = death_int_s(4)*r_m^3;
    
    if i == N_full + 1
        % bubbles leaving from free surface (m^3/(m^2 * s))
        
        % if we're looking at the bottom node, just take its value
        if i == 1
            death_term = 4/3 * pi * sum(r_q(i,:).^(3) .* w_q(i,:) .* (u_rise(i,:) - u_bulk) );
        else
            % above the bottom node, linearly interpolate/extrapolate to
            % get the value wherever the free surface is
            death_term_i = 4/3 * pi * sum(r_q(i,:).^(3) .* w_q(i,:) .* (u_rise(i,:) - u_bulk) );
            death_term_im1 = 4/3 * pi * sum(r_q(i-1,:).^(3) .* w_q(i-1,:) .* (u_rise(i-1,:) - u_bulk) );
            if node_level(i) < 0.5
                death_term = (0.5 + node_level(i))*death_term_i + (0.5 - node_level(i))*death_term_im1;
            else
                death_term_slope = (death_term_i - death_term_im1)/1;
                death_term = death_term_i + death_term_slope*( node_level(i) - 0.5 );
            end
        end
    else
        death_term = 0;
    end
    
    if i == 1
        % bubbles being convected out the bottom of the tank
        death_term_injector = 4/3 * pi * sum(r_q(i,:).^(3) .* w_q(i,:) .* (-u_rise(i,:) + u_bulk) );
        if death_term_injector < 0
            death_term_injector = 0;
        end
    end
    
    % mdot into bubbles from liquid
    %     mdot_bub_l(i) = V_l_star(i) * 4/3*pi * rho_tg_sat * (birth_term + growth_term);
    mdot_bub_l(i) = node_level(i) * V_node * 4/3 * pi * rho_tg_sat * (birth_term(i) + growth_term);
    
    if i == 1 || i == N_nodes
        mdot_bub_l(i) = mdot_bub_l(i)/2;
    end
    
    % mdot into bubbles from tg (really from bubbles into tg, but have to
    % keep sign convention for mdot)
    
    mdot_bub_tg(i) = - pi/4*D^2 * rho_tg_sat * death_term;
    
    % mdot into the bubbles
    mdot_bub(i) = mdot_bub_l(i) + mdot_bub_tg(i);
    
    if isinf(mdot_bub)
        disp('inf problem')
    end
    
end

mdot_bub_injector = pi/4*D^2 * rho_tg_sat * death_term_injector;

% need to deal with dw_dt and dg_dt for the nodes that have left the liquid
for i = N_full+2:N_nodes
    
    ind_node = (i-1)*N_ab + [1:N_ab];
    
    dw_dt(ind_node) = zeros(1,N_ab);
    dg_dt(ind_node) = zeros(1,N_ab);
    
end

% net rate of change of gas mass (mdot_bub_tg is negative)
mdot_tg = sum( -mdot_bub_tg );

% net rate of change of liquid mass
mdot_l = - sum(mdot_bub_l) - mdot_out_liq;

% HT from wall to liquid
% Qdot_lw = Qdot('lw',T_l,T_lw(1),rho_l,m_l,D);

Pr = P/P_cr;
h_20_lw = 1e3 * exp( 0.3092*log(Pr)^3 + 1.649*log(Pr)^2 + 3.641*log(Pr) + 5.272);
n_lw = 1 - 0.3*Pr^0.3;

q_20 = 20e3;

if T_lw(1) > T_s
    q_lw = ( ( (T_lw(1) - T_s)*h_20_lw )^(1/n_lw) / q_20 ) ^(n_lw/(1 - n_lw));
else
    q_lw = 0;
end

A_l = 4*V_l/D + pi/4*D^2;
Qdot_lw = C_qdot_lw * q_lw*A_l;

% Qdot_lw = 0;
% net HT into liquid
Qdot_l = Qdot_lw;

% HT into gas from wall
Qdot_gw = Qdot('gw',T_tg,T_gw(1),rho_tg,m_tg,D);
% Qdot_gw = 0;
% net HT into gas
Qdot_tg = Qdot_gw;

% not sure if this is correct... should it include a rhodot term?
% probably!!!
Vdot_bub = (sum( mdot_bub ) - mdot_out_vap)/ rho_tg_sat;

% this isn't actually Udot, but Udot without the P*Vdot term (hence the i)

Udot_li = - mdot_out_liq*h_l - sum( mdot_bub_l )*((h_tg_sat - h_l) + h_l) + Qdot_l - P*Vdot_bub;

Udot_tgi = Qdot_tg - sum( mdot_bub_tg )*(h_tg_sat) ;

% du_drho_tg = dh_drho_tg + P/rho_tg^2  - 1/rho_tg * dP_drho_tg;

du_drho_l = dh_drho_l + P/rho_l^2  - 1/rho_l * dP_drho_l;

Vdot_l = solve_for_Vdot(Udot_tgi, mdot_tg, m_tg, ...
    Udot_li, u_l, mdot_l, m_l, du_drho_l, Cv_l, dP_drho_l, V_l, ...
    dP_dT_l, V_tg, P, drho_dx_P_tg, drho_dP_x_tg, u_tg_v_sat, ...
    u_tg_l_sat, x_tg, ...
    du_dT_sat_tg_v, du_dT_sat_tg_l, dP_dT_tg_sat, ...
    rho_tg_v, u_tg, guesses, Vdot_bub);

% if Vdot_l == pi
%     disp('vdot error')
%     constants.error_detected = 1;
%     Vdot_l = guesses.Vdot_l;
%     keyboard
% end

Vdot_tg = - Vdot_l - Vdot_bub;

Udot_tg = Udot_tgi - P*Vdot_tg;

Udot_l = Udot_li - P*Vdot_l;

rhodot_l = mdot_l/V_l - m_l/V_l^2 * Vdot_l;

Tdot_l = ( ( Udot_l - u_l*mdot_l )/m_l - du_drho_l*rhodot_l )/Cv_l;

% mass of wall exposed to liquid
m_lw = tank_wall_mass(V_l,D,rho_w,t_w);

% mass of wall exposed to gas
m_gw = tank_wall_mass(V_tg,D,rho_w,t_w);

% HT from air to gas wall
Qdot_agw = Qdot('agw',T_air,T_gw(end),rho_tg,m_tg,D);

% HT from air to liquid wall
Qdot_alw = Qdot('alw',T_air,T_lw(end),rho_l,m_l,D);

% conduction from liquid wall to gas wall
L_tank = 4*V_tank/(pi*D^2);
L_wc = L_tank/2;
Qdot_wc = k_w*(T_lw - T_gw)*pi*D*t_w/L_wc;

% rate of change of mass of gass wall
mdot_gw = 4*Vdot_tg*t_w*rho_w/D;

% rate of change of temperature of gas wall
% Tdot_gw = (Qdot_agw - Qdot_gw + Qdot_wc + cv_w*mdot_gw*(T_lw - T_gw))/(m_gw*cv_w);

% rate of change of temperature of liquid wall
% Tdot_lw = (Qdot_alw - Qdot_lw - Qdot_wc)/(m_lw*cv_w);

Tdot_lw = wall_conduction(T_lw, Qdot_lw, Qdot_alw, constants);
Tdot_gw = wall_conduction(T_gw, Qdot_gw, Qdot_agw, constants);


% all derivatives
% 1 = m_tg
% 2 = U_tg
% 3 = T_gw
% 4 = m_l
% 5 = T_l
% 6 = T_lw
% 7:(N + 6) = weights
% N+7 : 2N+6 = weighted abscissas

dy = [mdot_tg;
    Udot_tg;
    Tdot_gw;
    mdot_l;
    Tdot_l;
    Tdot_lw;
    dw_dt';
    dg_dt'];

if ~sum(isreal(dy(:))) || sum(isnan(dy(:)))
    disp('complex or nans in derivatives')
    error_flag = 1;
end

% correct for negative stuff
if ~isempty(ind_negative) && sum( dy(ind_negative) < 0 ) > 0
    % we have negative dydt where y < 0!
    for i = 1:length(ind_negative)
        if dy(ind_negative(i)) < 0
            dy(ind_negative(i)) = 0;
        end
    end
end


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
%     debug_data.diff_eqns_error_flag = error_flag;
    
    varargout{2} = debug_data;
end



