function varargout = bubble_growth(varargin)
%% define input parameters

switch computer
    case {'GLNXA64'}
        % on server
        addpath('~/refprop/refprop');
        t_save = 60; % save interval in minutes
        plot_stuff = 0;
        save_stuff = 1;
        
    case {'MACI64'}
        % on laptop
        t_save = 10; % save interval in minutes
        plot_stuff = 1;
        save_stuff = 1;
end


constants.nuc_model = 'SJ';

N_mom = 6;
IC_moments = gamma_dist_moments( 5e-10, 1e-3, N_mom);

alpha_ic = 1e-9;
alpha_mom = 4/3 * pi * IC_moments(4);

IC_moments = IC_moments * alpha_ic/alpha_mom;

[r_ic, w_ic] = PD_method(IC_moments);

% N_mom = 6;

% r_ic = [3e-6 2e-3 7e-4];
% w_ic = [1 0.01  1e-3];

% r_ic = [2e-5 2e-4 5e-4];
% w_ic = [1 0.01  1e-3];

% r_ic = [3e-8 2e-8 7e-8];
% w_ic = [1 0.01  1e-3];

DQMOM_IC = [w_ic(:); r_ic(:).*w_ic(:)];

% DQMOM_IC = [1-6; 1e-6; 1e-6; 2e-1; 1e-4; 3e-9];

N_ab = N_mom/2;
constants.N_ab = N_ab; % number of abscissas
% (2*N = number of moments, going from 0 to 2N-1)

constants.phi = 1e-3;

clock_save = clock;

% close all

if nargin == 0
    constants.C_rdot = 3;
    constants.C_nuc_rate = 3;
    constants.n_nuc_freq = 3;
    constants.checking_gauss_error = 0;
    constants.C_death_rate = 1;
    constants.C_erf_death_rate = 10;
    % constants.n_death_rate = 2;
    % constants.alpha_lim = 0.05;
    constants.r_death = 0.5 * 0.0254;
    
    %     E = 3*3.2e2;
    %     A_inj = 0.8*2.217e-5;
    specified_case = 6;
    if ~exist('A_inj','var')
        switch specified_case
            case 1
                A_inj = 9.347e-5;
                E = 1.3e3;
            case 2
                A_inj = 2.277e-5;
                E = 2.4e2;
            case 3
                A_inj = 2.936e-5;
                E = 5.3e2;
            case 4
                A_inj = 2.637e-5;
                E = 1e3;
            case 5
                A_inj = 1.56e-7;
                E = 3.4e2;
            case 6
                A_inj = 6.977e-7;
                E = 7.5e2;
        end
    end
else
    A_inj = varargin{1};
    specified_case = varargin{2};
    E = 1;
    constants.C_rdot = varargin{3};
    constants.C_nuc_rate = varargin{4};
    constants.n_nuc_freq = varargin{5};
    constants.checking_gauss_error = 0;
    constants.C_death_rate = varargin{6};
    constants.n_death_rate = varargin{7};
    constants.alpha_lim = varargin{8};
end


if specified_case == 0
    % N2O test 11 from my data
    Ti = 280.4;           % [K] initial temperature
    fill_level = 0.87;        % [] initial fill_level ratio (by volume)
    E = 2.1e4;          % [] heat transfer multiplier
    V_tank = 1.80e-4;   % [m^3] tank volume
    L_tank = 0.356;     % [m] tank length
    Cd = 1;         % [] injector Cd
    A_inj = 7e-7;       % [m^2] injector area
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

% initialize program parameters
h = 1e-12;           % [s] initial time step
running = 1;        % [] switch, 1 = program running, 0 = program stopped
rel_tol = 1e-5;     % [] max relative error allowed in adaptive scheme
abs_tol = 1e3;     % [] max absolute error allowed in adaptive scheme
min_error = 1e-3;   % [] min error (relative to error_tol) before step size is increased
h_max = 1;       % [s] max allowable time step
h_min = 1e-16;      % [s] min allowable time step
t_end = 300;         % [s] end time (if LRO doesn't happen first)
LRO_tol = 5e-3;     % [s] tolerance for resolving the LRO point
dT_sup_tol = 1e-14;

ode_solver = 'DP'; % [] options: 'RKF' for runge-kutta-fehlberg
% 'euler' for 1st order euler
% 'RK4' for 4th order runge-kutta
% 'CK' for cash-karp
% 'DP' for dormand-prince

N_dim = 6 + N_mom;

% fsolve_options = optimset('display','off');

load PDT_table

% need this code for using qinterp2:
Tvec_table = PDT.T;
Pvec_table = PDT.P;

[Tgrid_table, Pgrid_table] = meshgrid(Tvec_table, Pvec_table);
PDT.T = Tgrid_table;
PDT.P = Pgrid_table;

% for differential equation solver:
% y(n+1) = y(n) + sum(i = 1 -> s) of b(i)*k(i)
% k(i) = h*f(stuff)
% h = step size
% f = the ode's  (dy/dt)
% f for k(1) = f( t(n), y(n) )
% f for k(2) = f( t(n) + c(2)*h , y(n) + a(2,1)*k(1) )
% f for k(3) = f( t(n) + c(3)*h , y(n) + a(3,1)*k(1) + a(3,2)*k(2) )
% ...
% the bs denotes b_star, and is used for error estimation


%% initialize things

switch ode_solver
    case 'RK4'
        % 4th order runge-kutta
        adaptive = 0;
        
        a = [0  0   0   0;
            .5  0   0   0;
            0   .5  0   0;
            0   0   0   1];
        
        c = [0; .5; .5; 1];
        
        b = [1/6; 1/3; 1/3; 1/6];
        
    case 'euler'
        % 1st order euler
        adaptive = 0;
        
        a = [0];
        
        b = 1;
        
        c = [0];
        
    case 'RKF'
        adaptive = 1;
        
        % 4th/5th order runge kutta fehlberg
        % a, b, bs (ie b star), and c comprise the butcher tableau
        % c is what's usually shown on the vertical axis
        % b is for the 5th order solution, bs for the 4th
        a = [0          0           0           0           0       0;
            .25         0           0           0           0       0;
            3/32        9/32        0           0           0       0;
            1932/2197   -7200/2197  7296/2197   0           0       0;
            439/216     -8          3680/513    -845/4104   0       0;
            -8/27       2           -3544/2565  1859/4104   -11/40  0];
        
        c = [0;         .25;        3/8;        12/13;      1;      .5];
        
        b = [16/135;    0;          6656/12825; 28561/56430;-9/50;  2/55];
        
        bs = [25/216;   0;          1408/2565;  2197/4104;  -1/5;   0];
        
    case 'CK'
        adaptive = 1;
        
        % cash-karp
        c = [0; 1/5; 3/10; 3/5; 1; 7/8];
        
        b = [37/378      0           250/621     125/594     0               512/1771]';
        
        bs = [2825/27648	0           18575/48384	13525/55296	277/14336       1/4]';
        
        a = [0          0       0           0               0           0;
            1/5        0       0           0               0           0;
            3/40       9/40    0           0               0           0;
            3/10       -9/10   6/5         0               0           0;
            -11/54      5/2     -70/27      35/27           0           0;
            1631/55296	175/512 575/13824	44275/110592	253/4096	0];
        
    case 'DP'
        adaptive = 1;
        
        % dormand-prince
        
        c = [0              1/5     3/10        4/5         8/9         1           1]';
        
        b = [35/384         0       500/1113	125/192	-2187/6784      11/84       0]';
        
        bs = [5179/57600	0       7571/16695	393/640	-92097/339200	187/2100	1/40]';
        
        
        a =[0           0           0           0           0       0       0;
            1/5         0           0           0           0       0       0;
            3/40        9/40        0           0           0       0       0;
            44/45       -56/15      32/9        0           0       0       0;
            19372/6561	-25360/2187	64448/6561	-212/729    0       0       0;
            9017/3168	-355/33     46732/5247	49/176	-5103/18656 0       0;
            35/384      0           500/1113	125/192	-2187/6784	11/84   0];
        
        
        
        
end

constants.error_detected = 0;

s = length(c); % number of stages in the scheme

k = zeros(N_dim,s);
ti = 0;
n = 1;              % [] counter

t = 0;



%
[rho_l, rho_tg, P] = refpropm('+-P','T',Ti,'Q',0.5,'N2O');
[u_tg] = refpropm('U', 'T', Ti ,'Q', 1, 'N2O');
P = P*1e3;

% P = 1e3*refpropm('P','T',Ti,'Q',0,'N2O');
%
% [rho_tg_l, rho_tg_v, ~, u_tg_v] = n2o_fits_for_getting_P(P);
% u_tg = u_tg_v;
% % rho_l = rho_tg_l;
% rho_tg = rho_tg_v;
%
% rho_l = qinterp2(PDT.T, PDT.P, PDT.D_liq, Ti, P/1e3);


T_sat = Ti;

for i = 1:2*N_ab
    mom(i) = sum( r_ic.^(i-1) .* w_ic );
end

% net volume of bubbles, per unit volume of liquid
V_bubi = 4/3*pi*mom(4);

V_l_star = fill_level*V_tank;
V_bub = V_l_star * V_bubi;


V_tg = (1 - fill_level)*(V_tank - V_bub);

V_l = V_l_star*(1 - V_bubi);

m_l = V_l*rho_l;
m_tg = V_tg*rho_tg;

U_tg = m_tg * u_tg;

T_l = Ti;
% T_tg = Ti;

y(:,1) = [m_tg; U_tg; Ti; m_l; Ti; Ti];
y((N_dim - N_mom + 1):(N_dim), 1) = DQMOM_IC(:);

% 1 = m_tg
% 2 = U_tg
% 3 = T_gw
% 4 = m_l
% 5 = T_l
% 6 = T_lw
% 7:(N + 6) = weights
% N+7 : 2N+6 = weighted abscissas

[K_b, N_A, h_planck] = universal_constants('boltzmann', 'avagadro', 'planck');

[P_cr, T_cr] = refpropm('PT', 'C', 0, '', 0, 'N2O');
P_cr = P_cr*1e3;

constants.E = E;
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

guesses.P = P;
guesses.rho_tg = rho_tg;
guesses.rho_l = rho_l;
guesses.Vdot_l = 0;

V_bub = 0;
alpha = 0;

dT_superheat = 0;

% I don't seem to be using this term anymore so just leave it
derivatives = zeros(5,1);

% initialize f
f = zeros(N_dim,1);

min_flag = 0;
peak_flag = 0;

%% ODE solver

% begin looping
while running == 1;
    %     figure(1)
    %     hold on
    %     plot(t,P,'bo')
    %     if ~stop.requested
    
    starti = max([n-3, 1]);
    
    Pdot = bdiff(P,starti,n,t,adaptive);
    rhodot_l = bdiff(rho_l,starti,n,t,adaptive);
    rhodot_tg = bdiff(rho_tg,starti,n,t,adaptive);
    Vdot_l(n+1) = V_tank*bdiff(V_l/V_tank,starti,n,t,adaptive);
    Vdot_tg(n+1) = V_tank*bdiff(V_tg/V_tank,starti,n,t,adaptive);
    
    
    % check for min
    if (t(n) > 0.05) && (Pdot) > 0
        if min_flag == 0;
            min_flag = 1;
            t_min = t(n);
            n_min = n;
        end
    end
    
    % check for peak
    if (min_flag == 1) && (peak_flag == 0)
        if (abs(Pdot) < 5e3) && (t(n) > 1.25*t_min)
            peak_flag = 1;
            %             running = 0;
            n_peak = n;
            t_peak = t(n);
        end
    end
    
    if peak_flag == 1
        if t(n) > t_peak*1.5
            %             running = 0;
        end
    end
    
    if t(n) > 2
        %         running = 0;
    end
    
    
    mdot_l = f(4);
    rhodot_l = mdot_l/V_l(n) - m_l(n)/V_l(n)^2*Vdot_l(n);
    
    mdot_tg = f(1);
    rhodot_tg = mdot_tg/V_tg(n) - m_tg(n)/V_tg(n)^2*Vdot_tg(n);
    
    if h > 5*h_min
        
        derivatives = 0.5*[Pdot; rhodot_l; rhodot_tg; Vdot_l(n+1); Vdot_tg(n+1)] + 0.5*derivatives;
        
        guesses.Vdot_l = 0.5*Vdot_l(n+1) + 0.5*guesses.Vdot_l;
        
    end
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
        
        fprintf(['t = %#4.4g, dt = %#4.4g, P = %#4.4g, V_bub = %#4.4g, dT_sup = %#6.6g,' ...
            'alpha = %#4.4g, T_l = %#4.4g, T_tg = %#4.4g, m_l = %#4.4g,'...
            'm_tg = %#4.4g, fill_level%% = %#4.4g, Vdot_l = %#4.4g, Vdot_tg = %#4.4g,'...
            ' rhodot_l = %#4.4g, rhodot_tg = %#4.4g\n'],...
            t(n), t(n) - t(max([1, n-1])), P(n)/6894.8, V_bub(n), dT_superheat(n), alpha(n), T_l(n), 0,...
            m_l(n), m_tg(n), 100*fill_level(n), Vdot_l(n+1),...
            Vdot_tg(n+1), rhodot_l, rhodot_tg);
        
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
        
        if (LRO_slope*h < -0.003/100) && (fill_level(n) < 0.1/100)
            h = h/4;
            
        elseif (LRO_slope*h < -0.03/100) && (fill_level(n) < 1/100);
            h = h/2;
        elseif (LRO_slope*h < -0.3/100) && (fill_level(n) < 5/100);
            h = h/2;
        end
        
        
        
        % also check if we're close to going subcooled -> superheated
        
        if dT_superheat(n) < 0
            
            % slope of dT_superheat curve
            dTs_slope = bdiff(dT_superheat,starti,n,t,adaptive);
            
            % projected t_sup (t when dT_sup = 0)
            t_sup = -dT_superheat(n)/dTs_slope + t(n);
            
            h_sup = t_sup - t(n); % distance to t_sup
            
            % if the step we're about to take is >3/4 the distance to LRO
            % and the distance te LRO is bigger than the tolerance
            if (h > 2*h_sup && h_sup > dT_sup_tol) && (h_sup > 0);
                
                % set h to 1/2 the distance to LRO (ie refine)
                h = 0.5*h_sup;
                
            end
            
            if (dTs_slope*h < -0.003/100) && (dT_superheat(n) > -0.1/100)
                h = h/4;
                
            elseif (dTs_slope*h < -0.03/100) && (dT_superheat(n) > -1/100);
                h = h/2;
            elseif (dTs_slope*h < -0.3/100) && (dT_superheat(n) > -5/100);
                h = h/2;
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
            
            % 1 = m_tg
            % 2 = U_tg
            % 3 = T_gw
            % 4 = m_l
            % 5 = T_l
            % 6 = T_lw
            % 7:(N + 6) = weights
            % N+7 : 2N+6 = weighted abscissas
            
            
            if i == 1
                
                % f = f( t(n) , y(n) )
                
                f = diffeqns(y(:,n), constants, guesses, PDT);
            else
                
                % f for k(2) = f( t(n) + c(2)*h , y(n) + a(2,1)*k(1) )
                % f for k(3) = f( t(n) + c(3)*h , y(n) + a(3,1)*k(1) + a(3,2)*k(2) )
                % and so on
                
                a_k_term = sum( (ones(N_dim,1)*a(i,1:i-1)).*k(:,1:i-1) ,2 );
                
                y_new = y(:,n) + a_k_term;
                
                %                 mom_part = ((N_dim - N_mom + 1):N_dim);
                %
                %                 if sum( y_new(mom_part) < 0) > 0
                % %                     disp('negative moments')
                %                 end
                
                f = diffeqns(y_new, ...
                    constants, guesses, PDT);
            end
            
            k(:,i) = f*h;
            
            %             y_new = (y(:,n) + sum( (ones(N_dim,1)*a(i,1:i)).*k(:,1:i),2 ) );
            
            %             % check for high T
            %             if (y_new(5) > T_cr)
            %                 disp('problem - temperature went above critical')
            %                 error_flag = 1;
            %                 y_new = y(:,n);
            %             end
            
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
        
        % 1 = m_tg
        % 2 = U_tg
        % 3 = T_gw
        % 4 = m_l
        % 5 = T_l
        % 6 = T_lw
        % 7:(N + 6) = weights
        % N+7 : 2N+6 = weighted abscissas
        
        
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
            
            [rel_err, ind_max_rel_err] = max(rel_err(isfinite(rel_err)));  % fix rel_err to the maximum finite value of rel_err
            
            
            
            abs_err = abs(err);
            
            %             abs_err = abs_err(1:6); % remove the bubble distribution terms
            
            abs_err = max(abs_err(isfinite(abs_err)));  % do the same for abs_err
            
            % check for possible problems: isempty statements are in case
            % abs and rel err are both full of non-finite values
            % isnan checks for nan's
            % isreal checks for imaginary numbers
            error_conditions = isempty(rel_err) + ...
                isempty(abs_err) +  ...
                isnan(sum(err)) + ...
                ~isreal(sum(y(:,n+1))) + ...
                (y(5,n+1) > T_cr) + ...
                error_flag;
            
            
            % 1 = m_tg
            % 2 = U_tg
            % 3 = T_gw
            % 4 = m_l
            % 5 = T_l
            % 6 = T_lw
            % 7:(N + 6) = weights
            % N+7 : 2N+6 = weighted abscissas
            
            % if any of those fail, set rel_err large so that the step gets
            % recomuputed
            if error_conditions > 0
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
    
    % 1 = m_tg
    % 2 = U_tg
    % 3 = T_gw
    % 4 = m_l
    % 5 = T_l
    % 6 = T_lw
    % 7:(N + 6) = weights
    % N+7 : 2N+6 = weighted abscissas
    
    m_tg(n+1) = y(1,n+1);
    U_tg(n+1) = y(2,n+1);
    m_l(n+1) = y(4,n+1);
    T_l(n+1) = y(5,n+1);
    
    %     disp('outside of loop')
    
    w_q = y(7:(N_ab+6),n+1);
    g_q = y(N_ab+7:(2*N_ab+6),n+1);
    
    r_q = g_q./w_q;
    
    for i = 1:2*N_ab
        mom(i,n+1) = sum( r_q.^(i-1) .* w_q );
    end
    
    
    if sum(mom<0) > 0
        disp('negative moments')
    end
    
    % net volume of bubbles, per unit volume of liquid
    V_bubi = 4/3*pi*mom(4,n+1);
    
    % get system pressure
    P(n+1) = get_P_from_mU_mT(m_tg(n+1), U_tg(n+1), m_l(n+1), T_l(n+1), ...
        V_tank, V_bubi, PDT, guesses);
    
    
    % saturation temp based on pressure
    T_sat(n+1) = refpropm('T','P',P(n+1)/1e3,'Q',0.5,'N2O');
    
    % liquid superheat
    dT_superheat(n+1) = T_l(n+1) - T_sat(n+1);
    
    
    % calculate liquid and vapor density based on T's and P
    %     rho_l(n+1) = get_D_from_TP(T_l(n+1), P(n+1), guesses);
    
    
    rho_l(n+1) = qinterp2(PDT.T,PDT.P,PDT.D_liq,T_l(n+1),P(n+1)/1e3);
    rho_tg(n+1) = refpropm('D', 'P', P(n+1)/1e3, 'Q', 1, 'N2O');
    
    guesses.rho_l = rho_l(n+1);
    
    
    % get volumes based on mass and density
    V_l(n+1) = m_l(n+1)/rho_l(n+1);
    V_tg(n+1) = m_tg(n+1)/rho_tg(n+1);
    
    V_l_star(n+1) = V_l(n+1)/(1 - V_bubi);
    
    
    % fill level based on liquid and tank volume
    fill_level(n+1) = V_l_star(n+1)/V_tank;
    
    % actual volume of all the bubbles
    V_bub(n+1) = V_bubi * V_l(n+1)/(1 - V_bubi);
    
    
    T_lw = y(6,n+1);
    
    Qdot_lw(n+1) = Qdot('lw',T_l(n+1), T_lw ,rho_l(n+1), m_l(n+1),D);
    
    %     V_bubi = 4/3*pi*y(7+3,n); % bubble volume per volume of liquid
    alpha(n+1) = V_bubi;
    
    
    guesses.P = P(n+1);
    guesses.rho_tg = rho_tg(n+1);
    guesses.rho_l = rho_l(n+1);
    
    t(n+1) = t(n) + h;
    
    n = n + 1;
    
    if (sum(abs(imag(y(:,n)))) > 0) || (sum(isnan(y(:,n))) > 0)
        running = 0;
        disp('imaginary or nans')
    end
    
    if (t(n) > t_end) || ( fill_level(n) < 1e-6)
        %         disp('reached end t or ran out of liquid')
        running = 0;
    end
    
    if t(n) > ti
        ti = ti + 0.1;
        %         disp(num2str(t(n)))
    end
    
    if P(end) < 2e5
        running = 0;
        disp('pressure got low')
    end
    
    alpha_next = alpha(end) + (alpha(end) - alpha(end-1));
    
    if alpha_next > 0.95
        running = 0;
        disp('alpha -> 1')
    end
    
    clock_now = clock;
    
    time_since_save = etime(clock_now, clock_save);
    
    if time_since_save >= t_save*60;
        if save_stuff == 1
            save('bubble_growth_sim_data','-v7.3')
        end
        clock_save = clock;
    end
    
    
    if constants.error_detected
        runnung = 0;
    end
    
    %         if rem(n,n_save) == 0
    %             save('bubble_growth_sim_data')
    %         end
    
    %         figure(1)
    %         subplot(1,2,1)
    %         hold on
    %         plot(t(n), P(n), 'ks')
    %
    %         subplot(1,2,2)
    %         loglog(r_dist, N_r + 1e-100)
    %         set(gca, 'xscale', 'log')
    %         set(gca, 'yscale', 'log')
    
    %     else
    %         stop.close
    %         error('you stopped me')
    %
    %     end
    
end

%% plotting and output
if nargout > 0
    %     P_LRO = P(end);
    %     T_LRO = T_l(end);
    %     t_LRO = t(end);
    %     varargout{1} = P_LRO/P(1);
    %     varargout{2} = T_LRO/T_l(1);
    %     varargout{3} = t_LRO;
    %
    %     if nargout > 3
    %         varargout{4} = P;
    %         varargout{5} = T_l;
    %         varargout{6} = t;
    %     end
    
    % if exist('n_min','var')
    %     varargout{1} = t(n_min);
    % varargout{3} = P(n_min);
    %
    % else
    %
    %     varargout{1} = NaN;
    %     varargout{3} = NaN;
    % end
    %
    % if exist('n_peak','var')
    %
    % varargout{2} = t(n_peak);
    % varargout{4} = P(n_peak);
    % else
    %
    % varargout{2} = NaN;
    % varargout{4} = NaN;
    % end
    
    varargout{1} = t;
    varargout{2} = P;
    
else
    
    if save_stuff == 1
        
        save('bubble_growth_sim_data','-v7.3')
        
    end
    
    % 1 = m_tg
    % 2 = U_tg
    % 3 = T_gw
    % 4 = m_l
    % 5 = T_l
    % 6 = T_lw
    % 7:(N + 6) = weights
    % N+7 : 2N+6 = weighted abscissas
    
    
    w_q = y(7:(N_ab+6),:);
    g_q = y(N_ab+7:(2*N_ab+6),:);
    
    r_q = g_q./w_q;
    for j = 1:n
        for i = 1:2*N_ab
            mom(i,j) = sum( r_q(:,j).^(i-1) .* w_q(:,j) );
        end
    end
    
    V_bubi = 4/3*pi*mom(4,:);
    A_bubi = 4*pi*mom(3,:);
    
    V_bub = V_bubi.*V_l_star;
    A_bub = A_bubi.*V_l_star;
    
    if plot_stuff == 1
        %
        figure(1)
        hold on
        plot(t,P/1e6,'k-')
        xlabel('Time [s]')
        ylabel('Pressure [MPa]')
        hold on
        if exist('t_peak','var')
            plot(t_min, P(n_min)/1e6, 'ko')
            plot(t_peak, P(n_peak)/1e6, 'ks')
        end
        
        figure(2)
        hold on
        plot(t,y(5,:),'k-',t,T_sat,'r:')
        legend('Liquid','T_{sat}(P)')
        ylabel('Temperature')
        xlabel('Time [s]')
        legend('Vapor','T_{sat}(P)')
        
        figure(3)
        hold on
        plot(t,y(6,:),'k-',t,y(3,:),'b--')
        title('wall temp')
        xlabel('Time [s]')
        legend('liquid','vapor')
        title('wall temp')
        
        figure(4)
        hold on
        plot(t,y(4,:),'k-',t,y(1,:),'b--')
        title('Mass')
        xlabel('Time [s]')
        legend('Liquid','Vapor')
        title('masses')
        
        figure(5)
        hold on
        plot(t,mom)
        xlabel('Time [s]')
        ylabel('Moments [various]')
        set(gca,'yscale','log')
        
        figure(7)
        hold on
        plot(t, A_bub,'k')
        xlabel('Time [s]')
        ylabel('A/V [1/m]')
        title('interfacial area per volume')
        
        figure(8)
        hold on
        plot(t, V_bubi,'k')
        xlabel('Time [s]')
        ylabel('gas holdup')
        
        figure(9)
        hold on
        plot(t, 6*V_bub./A_bub, 'k')
        xlabel('Time [s]')
        ylabel(' sauter diameter [m]')
        set(gca,'yscale','log')
        
        figure(10)
        hold on
        plot(t, r_q)
        xlabel('Time [s]')
        ylabel('abscissas [m]')
        set(gca,'yscale','log')
        
        figure(11)
        hold on
        plot(t, w_q)
        xlabel('Time [s]')
        ylabel('weights [?]')
        set(gca,'yscale','log')
        
        
    end
    
    
end

%% differential equations
function dy = diffeqns(y, constants, guesses, PDT)


% retrieve constants
E = constants.E;
D = constants.D;
t_w = constants.t_w;
rho_w = constants.rho_w;
cv_w = constants.cv_w;
Cd = constants.Cd;
A_inj = constants.A_inj;
Po = constants.Po;
T_air = constants.T_air;
V_tank = constants.V_tank;
h = constants.h;
k_w = constants.k_w;
K_b = constants.K_b;
N_A = constants.N_A;
P_cr = constants.P_cr;
T_cr = constants.T_cr;
h_planck = constants.h_planck;
nuc_model = constants.nuc_model;
g = constants.g;
C_hamaker = constants.C_hamaker;

n_nuc_freq = constants.n_nuc_freq;
% n_death_rate = constants.n_death_rate;
C_death_rate = constants.C_death_rate;
% alpha_lim = constants.alpha_lim;

phi = constants.phi;
C_rdot = constants.C_rdot;
C_nuc_rate = constants.C_nuc_rate;

C_erf_death_rate = constants.C_erf_death_rate;

r_death = constants.r_death;

N = constants.N_ab; % number of abscissas
% (2*N = number of moments, going from 0 to 2N-1)

% retrieve derivatives calculated with backwards differencing
% Pdot = derivatives(1);
% rhodot_l = derivatives(2);
% rhodot_tg = derivatives(3);
% Vdot_l = derivatives(4);
% Vdot_tg = derivatives(5);

% retrieve variables
% 1 = m_tg
% 2 = U_tg
% 3 = T_gw
% 4 = m_l
% 5 = T_l
% 6 = T_lw
% 7:(N + 6) = weights
% N+7 : 2N+6 = weighted abscissas

m_tg = y(1);
U_tg = y(2);
T_gw = y(3);
m_l = y(4);
T_l = y(5);
T_lw = y(6);
w_q = y(7:6+N);
g_q = y(7+N:6+2*N);

if isnan(sum(y)) || ~isreal(sum(y))
    disp('problem: nans or imaginary y')
end

r_q = g_q./w_q;

for i = 1:2*N
    mom(i) = sum( r_q.^(i-1) .* w_q );
end

% fprintf('abscissas = ')
% fprintf('%6.6g, ', r_q)
% fprintf('\n')
% fprintf('weights = ')
% fprintf('%6.6g, ', w_q)
% fprintf('\n')
% fprintf('L*w = ')
% fprintf('%6.6g, ', w_q.*r_q)
% fprintf('\n')

% bubble volume per unit volume of liquid/bubble mixture (hence the i)
V_bubi = 4/3*pi*mom(4);

% get system pressure
P = get_P_from_mU_mT(m_tg, U_tg, m_l, T_l, V_tank, V_bubi, PDT, guesses);

if P == pi
    disp('P error')
    constants.error_detected = 1;
    P = guesses.P;
end

% get density of liquid and ullage based on temperature and pressure

% rho_l = get_D_from_TP(T_l, P, guesses);

rho_l = qinterp2(PDT.T, PDT.P, PDT.D_liq, T_l, P/1e3);
% rho_tg = qqinterp2(PDT.T, PDT.P, PDT.D_vap, T_tg, P, 'linear');

% [rho_tg, T_tg] = refpropm('DT', 'P', P/1e3, 'U', U_tg/m_tg, 'N2O');

[rho_tg_l, rho_tg_v, u_tg_l, u_tg_v] = n2o_fits_for_getting_P(P);

u_tg = U_tg/m_tg;
x = (u_tg - u_tg_l)/(u_tg_v - u_tg_l);
alpha = 1/( 1 + rho_tg_v/rho_tg_l * (1 - x)/x );
rho_tg = alpha*rho_tg_v + (1 - alpha)*rho_tg_l;

V_l = m_l/rho_l;

V_l_star = V_l/(1 - V_bubi);

V_tg = m_tg/rho_tg;

% volume of all bubbles
% V_bub = V_l * V_bubi;

% liquid properties
[h_l, dh_drho_l, drho_dP_l, u_l, Cv_l, dP_dT_l, k_l, Cp_l, MW] = ...
    refpropm('H!RUO#LCM','T',T_l,'D&',rho_l,'N2O');
dP_drho_l = 1e3/drho_dP_l;
dP_dT_l = dP_dT_l*1e3;
alpha_l = k_l/(rho_l * Cp_l);


% a lot of these properties aren't needed until later, but by calculating
% them here I can reduce the number of refprop calls. Note that the next
% two refprop calls are for ullage properties (liquid and vapor parts)

% properties needed for Vdot calculation (liquid)
[u_tg_l_sat, rho_tg_l, dP_dT_tg_sat, drho_dP_T, drho_dT_P, dh_dT_P, dh_dP_T...
    ,sigma, h_l_sat, T_tg] = refpropm('UDERW(*IHT', 'P', P/1e3, 'Q', 0, 'N2O');
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
    T_s, h_tg_sat] = ...
    refpropm('UDRW(*TH', 'P', P/1e3, 'Q', 1, 'N2O');
drho_dP_T = drho_dP_T * 1e-3;
dh_dP_T = dh_dP_T * 1e-3;

du_dT_P = dh_dT_P + P/rho_tg_v^2 * drho_dT_P;
du_dP_T = dh_dP_T + 1/rho_tg_v + P/rho_tg_v^2 * drho_dP_T;

du_dT_sat_tg_v = du_dT_P + du_dP_T * dP_dT_tg_sat;
drho_dT_v_sat = drho_dT_P + drho_dP_T * dP_dT_tg_sat;
% drho_dP_sat_tg = drho_dP_T + drho_dT_P / dP_dT_sat;

drho_dx_P_tg = -rho_tg^2 *(1/rho_tg_v - 1/rho_tg_l);
drho_dP_x_tg = (1/dP_dT_tg_sat) * rho_tg^2 * ( x/rho_tg_v^2 * drho_dT_v_sat + ...
    (1-x)/rho_tg_l^2 * drho_dT_l_sat );


rho_tg_sat = rho_tg_v;
% % temp of saturated surface based on pressure (and h of sat. vapor)
% [T_s, h_tg_sat, rho_tg_sat] = refpropm('THD','P',P/1e3,'Q',1,'N2O');
%
% % saturated liquid enthalpy at P
% [sigma, h_l_sat] = refpropm('IH','P',P/1e3,'Q',0,'N2O');

% heat of vaporization (at saturation)
h_lv = h_tg_sat - h_l_sat;

% bubble calculations

% superheat = T - T_sat
deltaT_sup = T_l - T_s;

[P_sat, s_liq_sat, h_liq_sat] = refpropm('PSH', 'T', T_l, 'Q', 0, 'N2O');
P_sat = 1e3*P_sat;

% if superheated, then calculate bubble stuff
if deltaT_sup > 1e-6
    
    if sum(abs(imag([r_q(:); w_q(:)]))) > 0
        fprintf('imaginary abscissas or weights. moments:')
        fprintf('%0.6g\t',mom)
        fprintf('\n')
        
    end
    
    % jakob number
    Ja_T = Cp_l * rho_l * deltaT_sup/(rho_tg * h_lv);
    
    % bubble radius rate of change
    rdot = C_rdot * Ja_T^2 * alpha_l ./ r_q;
    
    % radius of new bubbles
    r_nuc = 2*sigma*T_s/(rho_tg * h_lv * deltaT_sup);
    
    % length of liquid volume [m]
    L_l = V_l_star / (pi * 0.25 * D^2);
    
    % surface area [m^2]
    A_l = pi * D * L_l + pi * 0.25 * D^2;
    
    switch constants.nuc_model
        
        case 'SJ'
            
            % rate of nucleation (mostly from shin & jones, 1993)
            
            % departure diameter [m]
            r_dep = 0.5 * 2.97e4 * (P/P_cr)^-1.09 * ( K_b * T_cr / (P_cr * MW) )^(1/3);
            % correlation from Jensen & Memmel, 1986
            % gives values on the order of 10-20 microns
            
            % non-dimensional cavity size (ie bubble nucleation size)
            r_c_star = r_nuc / r_dep;
            
            % non-dimensional nucleation rate
            N_ns_star = 1e-7 * r_c_star;
            
            % nucleation site density [1/m^2]
            nuc_density = N_ns_star * ( 0.5 / r_dep )^2;
            
            % nucleation frequency [Hz]
            nuc_freq = 1e4 * deltaT_sup^n_nuc_freq;
            
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
    
    
    
else
    
    % no superheat -> no bubble growth or nucleation
    r_nuc = 0;
    rdot = 0;
    nuc_rate = 0;
    
end

% generat scaled abscissas for more accurate calculations
r_m = max(abs(r_q));    % value to scale by

r_s = r_q/r_m;  % scaled abscissas

growth_int_s(1) = 0;

for i = 2:N*2
    growth_int_s(i) = (i-1)*sum( (1/r_m) * r_s.^(i-2) .* w_q .* rdot );
end

spec_nuc_rate = nuc_rate / V_l_star;
% nucleation rate per volume

% birth due to nucleation
% death due to bubbles getting too large
for i = 1:N*2
    birth_int_s(i) = (r_nuc/r_m).^(i-1) * spec_nuc_rate;
    death_int_s(i) = sum( r_s.^(i-1) .* w_q .* C_death_rate ...
        .* 0.5.*(1 + erf( C_erf_death_rate * (abs(r_q)/r_death - 1) ) ) );
end

% birth and death due to coalescence

nu_t = 0.0536*D^1.77/rho_l;
U_max = ( (1 - 0.75*V_bubi)/(1 - V_bubi) )*V_bubi * D^2/(48*nu_t);
mean_shear = 0.53*U_max/(0.5*D);

P2 = P;

liquid_height = V_l/(pi*D^2/4);

P1 = P2 + rho_l*g*liquid_height;

L = V_tank/(pi*D^2/4);

Q = L/8;

turb_diss = Q*g * P2*log(P1/P2) / ( pi * (0.5*D)^2 * (P1 - P2) );

coal_birth_s = zeros(2*N,1);
coal_death_s = zeros(2*N,1);

for k = 1:N*2
    for i = 1:N
        for j = 1:N
            rbi = r_q(i);
            rbj = r_q(j);
            dbi = 2*rbi;
            dbj = 2*rbj;
            
            qLS = 4/3*(rbi + rbj)^3*mean_shear;
            
            u_ri = sqrt( (2.14*sigma/(rho_l*dbi)) + 0.505*g*dbi);
            u_rj = sqrt( (2.14*sigma/(rho_l*dbj)) + 0.505*g*dbj);
            Sij = pi/4*(rbi + rbj)^2;
            qB = Sij*abs(u_ri - u_rj);
            
            qT = 0.089*pi*(dbi + dbj)^2 * turb_diss^(1/3) * sqrt(dbi^(2/3) + dbj^(2/3));
            
            % coalescence kernel
            
            rb_eq = ( (1/rbi + 1/rbj)/2 )^-1;
            
            t_cont = rb_eq^(2/3) / turb_diss^(1/3);
            
            film_i = 1e-4;
            film_f = (C_hamaker * rb_eq/(8*pi*sigma))^(1/3);
            t_coal = sqrt( rb_eq^3 * rho_l/(16 * sigma) ) * log( film_i / film_f);
            
            beta = (qT + qB + qLS)*exp( - t_coal / t_cont);
            
            coal_birth_s(k) = coal_birth_s(k) + 0.5 * w_q(i) * w_q(j) * ( r_s(i)^3 + r_s(j)^3 )^((k-1)/3) * beta;
            coal_death_s(k) = coal_death_s(k) + r_s(i)^(k-1) * w_q(i) * w_q(j) * beta;
        end
    end
end



    
    
dmom_dt_s = birth_int_s(:) - death_int_s(:) + growth_int_s(:) + coal_birth_s(:) - coal_death_s(:);


% if deltaT_sup > 0
%     fprintf('coalescence birth:\n')
%     fprintf('%6.6g\t', coal_birth_s./dmom_dt_s)
%         fprintf('\ncoalescence death:\n')
%     fprintf('%6.6g\t', coal_death_s./dmom_dt_s)
%     fprintf('\n')
% end

ill_conditioned = 1;

n = 1;
while ill_conditioned
    
    if n > 1
        
        r_sp = r_s + (0.5 - rand(size(r_s))).*abs(r_s)*1e-6;
        
    else
        r_sp = r_s;
    end
    
    for j = 0:(2*N - 1)
        if j == 0
            A1(1,:) = ones(1,N);
            A2(1,:) = zeros(1,N);
        elseif j == 1
            A1(2,:) = zeros(1,N);
            A2(2,:) = ones(1,N);
        else
            A1(j+1,:) = (1 - j)*r_sp.^j;
            A2(j+1,:) = j*r_sp.^(j);
        end
    end
    
    A = [A1 A2];
    
    cond_val = rcond(A);
    
    if (cond_val < 1e-15) || isnan(cond_val)
        ill_conditioned = 1;
    else
        ill_conditioned = 0;
        if n > 1
            %             disp('fixed')
        end
    end
    
    if n > 100
        ill_conditioned = 0;
        disp('gave up')
    end
    
    n = n + 1;
end

beta_q = dmom_dt_s;

alpha_q = A\beta_q;

a_q = alpha_q(1:N);
b_q_s = alpha_q(N+1:end);
b_q = b_q_s.*r_q;

dw_dt = a_q;
dg_dt = b_q;

birth_term = birth_int_s(4)*r_m^3;
growth_term = growth_int_s(4)*r_m^3;
death_term = death_int_s(4)*r_m^3;

% mdot into bubbles from liquid
mdot_bub_l = V_l_star * 4/3*pi * rho_tg_sat * (birth_term + growth_term);

% mdot into bubbles from tg
mdot_bub_tg = - V_l_star * 4/3*pi * rho_tg_sat * death_term;

mdot_bub = mdot_bub_l + mdot_bub_tg;

if isinf(mdot_bub)
    disp('inf problem')
end

% mass flow rate out via injector
mdot_out = A_inj*Cd*dyer_flow(Po, P, T_l, rho_l, P_sat, s_liq_sat, h_liq_sat);

% net rate of change of gas mass
mdot_tg = -mdot_bub_tg;

% net rate of change of liquid mass
mdot_l = - mdot_bub_l - mdot_out;

% HT from wall to liquid
Qdot_lw = Qdot('lw',T_l,T_lw,rho_l,m_l,D);

% net HT into liquid
Qdot_l = Qdot_lw;

% HT into gas from wall
Qdot_gw = Qdot('gw',T_tg,T_gw,rho_tg,m_tg,D);

% net HT into gas
Qdot_tg = Qdot_gw;

% this isn't actually Udot, but Udot without the P*Vdot term (hence the i)

Udot_li = - mdot_out*h_l - mdot_bub_l*(h_lv + (h_l_sat - h_l) + h_l) + Qdot_l;

Udot_tgi = Qdot_tg;

% du_drho_tg = dh_drho_tg + P/rho_tg^2  - 1/rho_tg * dP_drho_tg;

du_drho_l = dh_drho_l + P/rho_l^2  - 1/rho_l * dP_drho_l;

% not sure if this is correct... should it include a rhodot term?
Vdot_bub = mdot_bub / rho_tg_sat;

Vdot_l = solve_for_Vdot(Udot_tgi, mdot_tg, m_tg, ...
    Udot_li, u_l, mdot_l, m_l, du_drho_l, Cv_l, dP_drho_l, V_l, ...
    dP_dT_l, V_tg, P, drho_dx_P_tg, drho_dP_x_tg, u_tg_v_sat, ...
    u_tg_l_sat, x, ...
    du_dT_sat_tg_v, du_dT_sat_tg_l, dP_dT_tg_sat, guesses, Vdot_bub);

if Vdot_l == pi
    disp('vdot error')
    constants.error_detected = 1;
    Vdot_l = guesses.Vdot_l;
end


Vdot_tg = - Vdot_l - Vdot_bub;

Udot_tg = Udot_tgi - P*Vdot_tg;

Udot_l = Udot_li - P*Vdot_l;

% rhodot_tg = mdot_tg/V_tg - m_tg/V_tg^2 * Vdot_tg;

rhodot_l = mdot_l/V_l - m_l/V_l^2 * Vdot_l;

% Tdot_tg = ( ( Udot_tg - u_tg*mdot_tg )/m_tg - du_drho_tg*rhodot_tg )/Cv_tg;

Tdot_l = ( ( Udot_l - u_l*mdot_l )/m_l - du_drho_l*rhodot_l )/Cv_l;

% mass of wall exposed to liquid
m_lw = tank_wall_mass(V_l,D,rho_w,t_w);

% mass of wall exposed to gas
m_gw = tank_wall_mass(V_tg,D,rho_w,t_w);

% HT from air to gas wall
Qdot_agw = Qdot('agw',T_air,T_gw,rho_tg,m_tg,D);

% HT from air to liquid wall
Qdot_alw = Qdot('alw',T_air,T_lw,rho_l,m_l,D);

% conduction from liquid wall to gas wall
L_tank = 4*V_tank/(pi*D^2);
L_wc = L_tank/2;
Qdot_wc = k_w*(T_lw - T_gw)*pi*D*t_w/L_wc;

% rate of change of mass of gass wall
mdot_gw = 4*Vdot_tg*t_w*rho_w/D;

% rate of change of temperature of gas wall
Tdot_gw = (Qdot_agw - Qdot_gw + Qdot_wc + cv_w*mdot_gw*(T_lw - T_gw))/(m_gw*cv_w);

% rate of change of temperature of liquid wall
Tdot_lw = (Qdot_alw - Qdot_lw - Qdot_wc)/(m_lw*cv_w);

% all derivatives
% 1 = m_tg
% 2 = U_tg
% 3 = T_gw
% 4 = m_l
% 5 = T_l
% 6 = T_lw
% 7:(N + 6) = weights
% N+7 : 2N+6 = weighted abscissas
%
% figure(1)
% bar(r_q,w_q,0.1)

dy = [mdot_tg;
    Udot_tg;
    Tdot_gw;
    mdot_l;
    Tdot_l;
    Tdot_lw;
    dw_dt(:);
    dg_dt(:)];





