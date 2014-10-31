function varargout = ZK(varargin)

if nargin == 0
%     E = 3*3.2e2;
%     A_inj = 0.8*2.217e-5;
    specified_case = 2;
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
    E = varargin{3};
end

% A_inj
% tic
% problem parameter
% specified_case = 6;
% A_inj = 3.2e-6;
% plotting_LRO = 1;

% specified_case = 4; % [] 0 = no specified case, use following
% 1 = greg's setup from his JPC paper on modeling
% 2 = my tank testing experimental setup

% initialize program parameters
h = 1e-12;           % [s] initial time step
running = 1;        % [] switch, 1 = program running, 0 = program stopped
rel_tol = 1e-5;     % [] max relative error allowed in adaptive scheme
abs_tol = 1e3;     % [] max absolute error allowed in adaptive scheme
min_error = 1e-2;   % [] min error (relative to error_tol) before step size is increased
adaptive = 1;       % [] switch, 1 = use 4/5 Runge-Kutta_Fehlberg, 0 = plain 4th Runge-Kutta
h_max = 1e-0;       % [s] max allowable time step
h_min = 1e-16;      % [s] min allowable time step
t_end = 300;         % [s] end time (if LRO doesn't happen first)
LRO_tol = 5e-3;     % [s] tolerance for resolving the LRO point
non_adaptive_scheme = 1; % [] switch, 1 = 4th order Runge-Kutta, 0 = 1st order Euler

N_dim = 6;

fsolve_options = optimset('display','off');

load PDT_table

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

if specified_case == 0
            % can play with this one
        Ti = 292.5;           % [K] initial temperature
        fill_level = 0.9;        % [] initial fill_level ratio (by volume)
        %         E = 2e2;          % [] heat transfer multiplier
        V_tank = 0.011276;   % [m^3] tank volume
        L_tank = 62.3*0.0254;     % [m] tank length
        Cd = 0.5;         % [] injector Cd
        %         A_inj = 3.15e-7;       % [m^2] injector area
        Po = 1e5;           % [Pa] external pressure
        T_air = 293;        % [K] air temperature
        rho_w = 2700;       % [kg/m^3] density of wall material (aluminum)
        cv_w = 904;         % [J/kg.K] specific heat of wall (aluminum)
        t_w = 0.0254*1/8;   % [m] wall thickness
        D = sqrt(4/pi*V_tank/L_tank);
        % [m] tank diameter
        k_w = 167;          % [W/m.K] thermal conductivity of wall
else

[Ti, fill_level, V_tank, L_tank, ...
    Cd, Po, T_air, rho_w, cv_w, t_w, D, k_w] = initial_conditions(specified_case);

end
if adaptive == 0
    
    if non_adaptive_scheme
        % 4th order RK:
        a = [0  0   0   0;
            .5  0   0   0;
            0   .5  0   0;
            0   0   0   1];
        
        c = [0; .5; .5; 1];
        
        b = [1/6; 1/3; 1/3; 1/6];
    else
        
        % 1st order Euler
        
        a = [0];
        
        b = 1;
        
        c = [0];
    end
    
else
    % 4th/5th order runge kutta fehlberg
    a = [0          0           0           0           0       0;
        .25         0           0           0           0       0;
        3/32        9/32        0           0           0       0;
        1932/2197   -7200/2197  7296/2197   0           0       0;
        439/216     -8          3680/513    -845/4104   0       0;
        -8/27       2           -3544/2565  1859/4104   -11/40  0];
    
    c = [0;     .25; 3/8;           12/13;          1;      .5];
    
    b = [16/135; 0;  6656/12825;    28561/56430;    -9/50;  2/55];
    
    bs = [25/216; 0; 1408/2565;     2197/4104;      -1/5;   0];
    
end

s = length(c);

k = zeros(N_dim,s);
ti = 0;
n = 1;              % [] counter

t = 0;

[rho_l, rho_tg, P] = refpropmJEZ('+-P','T',Ti,'Q',0.5,'N2O');
P = P*1e3;

T_sat = Ti;

V_l = fill_level*V_tank;
V_tg = (1 - fill_level)*V_tank;

m_l = fill_level*V_tank*rho_l;
m_tg = (1 - fill_level)*V_tank*rho_tg;

T_l = Ti;
T_tg = Ti;

y(:,1) = [m_tg; Ti; Ti; m_l; Ti; Ti];

% 1 = m_tg
% 2 = T_tg
% 3 = T_gw
% 4 = m_l
% 5 = T_l
% 6 = T_lw

constants = [E; D; t_w; rho_w;
    cv_w; Cd; A_inj; Po;
    T_air; V_tank; 0; k_w];

guesses = [P; rho_tg; rho_l; 0];
P_guess = P;

derivatives = zeros(5,1);

f = zeros(6,1);

% begin looping
while running == 1;
    
    starti = max([n-3, 1]);
    
    Pdot = bdiff(P,starti,n,t,adaptive);
    rhodot_l = bdiff(rho_l,starti,n,t,adaptive);
    rhodot_tg = bdiff(rho_tg,starti,n,t,adaptive);
    Vdot_l(n+1) = V_tank*bdiff(V_l/V_tank,starti,n,t,adaptive);
    Vdot_tg(n+1) = V_tank*bdiff(V_tg/V_tank,starti,n,t,adaptive);
    
    mdot_l = f(4);
    rhodot_l = mdot_l/V_l(n) - m_l(n)/V_l(n)^2*Vdot_l(n);
    
    mdot_tg = f(1);
    rhodot_tg = mdot_tg/V_tg(n) - m_tg(n)/V_tg(n)^2*Vdot_tg(n);
    
    if h > 5*h_min
        
        derivatives = 0.5*[Pdot; rhodot_l; rhodot_tg; Vdot_l(n+1); Vdot_tg(n+1)] + 0.5*derivatives;
        
        guesses(4) = 0.5*Vdot_l(n+1) + 0.5*guesses(4);
        
    end
    %     derivatives = zeros(5,1);
    
    T_lw = y(6,n);
    
    Qdot_lw(n) = Qdot('lw',T_l(n), T_lw ,rho_l(n), m_l(n),D);
    
    if nargin == 0
        
        fprintf(['t = %8.6g, P = %6.4g, T_l = %6.4g, T_tg = %6.4g, m_l = %6.4g,'...
            'm_tg = %6.4g, fill_level%% = %6.4g, Vdot_l = %6.4g, Vdot_tg = %6.4g,'...
            'u_l = %6.4g, rhodot_l = %6.4g, rhodot_tg = %6.4g\n'],...
            t(n), P(n)/6894.8, T_l(n), T_tg(n), m_l(n), m_tg(n), 100*fill_level(n), Vdot_l(n+1),...
            Vdot_tg(n+1), y(5,end)/y(4,end),rhodot_l, rhodot_tg);
        
    end
    
    % if we're not on the first step and error is plenty small, increase
    % step size
    if  n > 1 && adaptive == 1
        
        % if error is < min_error
        if max(abs_err/abs_tol,rel_err/rel_tol) < min_error
            
            % make h bigger
            h = min(4*h,h_max);
            
        end
        
        % slope of fill level curve
        slope = bdiff(fill_level,starti,n,t,adaptive);
        %         slope = (fill_level(n) - fill_level(n-1))/(t(n) - t(n-1));
        
        % projected t_LRO
        t_LRO = -fill_level(n)/slope + t(n);
        
        h_LRO = t_LRO - t(n); % distance to t_LRO
        
        % if the step we're about to take is >3/4 the distance to LRO
        % and the distance te LRO is bigger than the tolerance
        if (h > 2*h_LRO && h_LRO > LRO_tol) && (h_LRO > 0);
            
            % set h to 1/2 the distance to LRO (ie refine)
            h = 0.5*h_LRO;
            
        end
        
        if (slope*h < -0.003/100) && (fill_level(n) < 0.1/100)
            h = h/4;
        
        elseif (slope*h < -0.03/100) && (fill_level(n) < 1/100);
            h = h/2;
        elseif (slope*h < -0.3/100) && (fill_level(n) < 5/100);
            h = h/2;
        end
        
        
        
        
        
    end
    
    error_OK = 0;
    
    
    while error_OK == 0
        % solving differential equations
        % i = counter for
        
        constants(11) = h;
        
        error_flag = 0;
        
        for i = 1:s
            
            % 1 = m_tg
            % 2 = T_tg
            % 3 = T_gw
            % 4 = m_l
            % 5 = T_l
            % 6 = T_lw
            
            
            if i == 1
                
                % f = f( t(n) , y(n) )
                
                f = diffeqns(y(:,n), constants, derivatives, guesses, PDT);
            else
                
                % f for k(2) = f( t(n) + c(2)*h , y(n) + a(2,1)*k(1) )
                % f for k(3) = f( t(n) + c(3)*h , y(n) + a(3,1)*k(1) + a(3,2)*k(2) )
                % and so on
                
                f = diffeqns(y_new, ...
                    constants, derivatives, guesses, PDT);
            end
            
            k(:,i) = f*h;
            
            y_new = (y(:,n) + sum( (ones(N_dim,1)*a(i,1:i)).*k(:,1:i),2 ) );
            
            if (y_new(2) > 305) || (y_new(5) > 305)
                disp('problem')
                error_flag = 1;
                y_new = y(:,n);
            end
            
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
        
        if adaptive == 1
            % using adaptive scheme, need to check error
            
            err = k*(b - bs);   % absolute error (diff. between 5th and 4th order estimates)
            
            rel_err = abs(err./( y(:,n) + 1e-6));  % relative error
            
            rel_err = max(rel_err(isfinite(rel_err)));  % fix rel_err to the maximum finite value of rel_err
            
            abs_err = max(abs(err(isfinite(err)))); % do the same for abs_err
            
            % check for possible problems: isempty statements are in case
            % abs and rel err are both full of non-finite values
            % isnan checks for nan's
            % isreal checks for imaginary numbers
            error_conditions = isempty(rel_err) + ...
                isempty(abs_err) +  ...
                isnan(sum(err)) + ...
                ~isreal(sum(y(:,n+1))) + ...
                (y(2,n+1) > 308) + ...
                (y(5,n+1) > 308) + ...
                error_flag;
            
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
                
                % not meeting error requirements
                
                if rel_err == 0 || abs_err == 0
                    % something odd happened, so reduce step size
                    
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
    % 2 = T_tg
    % 3 = T_gw
    % 4 = m_l
    % 5 = T_l
    % 6 = T_lw
    
    m_tg(n+1) = y(1,n+1);
    T_tg(n+1) = y(2,n+1);
    m_l(n+1) = y(4,n+1);
    T_l(n+1) = y(5,n+1);
    
    P(n+1) = get_P_from_mm_TT(m_tg(n+1), T_tg(n+1), m_l(n+1), T_l(n+1), V_tank, PDT, guesses);
    
    T_sat(n+1) = refpropm('T','P',P(n+1)/1e3,'Q',0.5,'N2O');
    
    % P_guess = guesses(1);
    rho_tg_guess = guesses(2);
    rho_l_guess = guesses(3);
    
    %     rho_l(n+1) = easy_D(P(n+1), T_l, PDT, rho_l_guess, fsolve_options);
    %     rho_tg(n+1) = easy_D(P(n+1), T_tg, PDT, rho_tg_guess, fsolve_options);
    %
    rho_l(n+1) = interp2(PDT.T,PDT.P,PDT.D_liq,T_l(n+1),P(n+1)/1e3,'linear');
    rho_tg(n+1) = interp2(PDT.T,PDT.P,PDT.D_vap,T_tg(n+1),P(n+1)/1e3,'linear');
    
    V_l(n+1) = m_l(n+1)/rho_l(n+1);
    V_tg(n+1) = m_tg(n+1)/rho_tg(n+1);
    
    fill_level(n+1) = V_l(n+1)/V_tank;
    
    guesses(1:3) = [P(n+1); rho_tg(n+1); rho_l(n+1)];
    
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
    
    
end

% toc
%
if nargout > 0
    P_LRO = P(end);
    T_LRO = T_l(end);
    t_LRO = t(end);
    varargout{1} = P_LRO/P(1);
    varargout{2} = T_LRO/T_l(1);
    varargout{3} = t_LRO;
    
        if nargout > 3
        varargout{4} = P;
        varargout{5} = T_l;
        varargout{6} = t;
    end
else
    
    figure(1)
    hold on
    plot(t,P/1e6,'k-')
    xlabel('Time [s]')
    ylabel('Pressure [s]')
    
    figure(2)
    hold on
    plot(t,y(5,:),'k-',t,y(2,:),'b--',t,T_sat,'r:')
    ylabel('Temperature')
    xlabel('Time [s]')
    
    legend('Liquid','Vapor','T_{sat}(P)')
    
    figure(3)
    hold on
    plot(t,y(6,:),'k-',t,y(3,:),'b--')
    title('wall temp')
    xlabel('Time [s]')
    
    legend('liquid','vapor')
    
    figure(4)
    hold on
    plot(t,y(4,:),'k-',t,y(1,:),'b--')
    title('Mass')
    xlabel('Time [s]')
    
    legend('Liquid','Vapor')
    
    figure(5)
    subplot(1,3,1)
    hold on
    plot(t,P/1e6,'k-')
    xlabel('Time [s]')
    ylabel('Pressure [MPa]')
    
    subplot(1,3,2)
    hold on
    plot(t,T_l,'k',t,T_tg,'b--')
    xlabel('Time [s]')
    ylabel('Temperature [K]')
    legend('Liquid','Vapor')
    
    subplot(1,3,3)
    hold on
    plot(t,m_l,'k-',t,m_tg,'b--')
    xlabel('Time [s]')
    ylabel('Mass [kg]')
    legend('Liquid','Vapor')
    
    figure(6)
    hold on
    plot(t(1:end-1),Qdot_lw./m_l(1:end-1),'k-')
    xlabel('Time [s]')
    ylabel('Qdot/m')
    
    
end

function dy = diffeqns(y, constants, derivatives, guesses, PDT)

% retrieve constants
E = constants(1);
D = constants(2);
t_w = constants(3);
rho_w = constants(4);
cv_w = constants(5);
Cd = constants(6);
A_inj = constants(7);
Po = constants(8);
T_air = constants(9);
V_tank = constants(10);
h = constants(11);
k_w = constants(12);

% retrieve derivatives calculated with backwards differencing
% Pdot = derivatives(1);
% rhodot_l = derivatives(2);
% rhodot_tg = derivatives(3);
% Vdot_l = derivatives(4);
% Vdot_tg = derivatives(5);

% retrieve variables
% 1 = m_tg
% 2 = T_tg
% 3 = T_gw
% 4 = m_l
% 5 = T_l
% 6 = T_lw
m_tg = y(1);
T_tg = y(2);
T_gw = y(3);
m_l = y(4);
T_l = y(5);
T_lw = y(6);
% 
% P_guess = guesses(1);
% rho_tg_guess = guesses(2);
% rho_l_guess = guesses(3);
Vdot_guess = guesses(4);

if isnan(sum(y)) || ~isreal(sum(y))
    disp('problem')
end

if T_tg > 305
    disp('problem')
end

P = get_P_from_mm_TT(m_tg, T_tg, m_l, T_l, V_tank, PDT, guesses);

% rho_l = easy_D(P, T_l, PDT, rho_l_guess, fsolve_options);

% rho_tg = easy_D(P, T_tg, PDT, rho_tg_guess, fsolve_options);

rho_l = interp2(PDT.T, PDT.P, PDT.D_liq, T_l, P/1e3, 'linear');
rho_tg = interp2(PDT.T, PDT.P, PDT.D_vap, T_tg, P/1e3, 'linear');

V_l = m_l/rho_l;

V_tg = m_tg/rho_tg;

% liquid properties
[h_l, dh_drho_l, drho_dP_l, u_l, Cv_l, dP_dT_l] = ...
    refpropmJEZ('H!RUO#','T',T_l,'D&',rho_l,'N2O');
dP_drho_l = 1e3/drho_dP_l;
dP_dT_l = dP_dT_l*1e3;

% gas properties
[h_tg, dh_drho_tg, drho_dP_tg, u_tg, Cv_tg, dP_dT_tg] = ...
    refpropmJEZ('H!RUO#','T',T_tg,'D&',rho_tg,'N2O');
dP_drho_tg = 1e3/drho_dP_tg;
dP_dT_tg = dP_dT_tg*1e3;

% temp of saturated surface based on pressure (and h of sat. vapor)

[T_s, h_tg_sat] = refpropmJEZ('TH','P',P/1e3,'Q',1,'N2O');

[P_sat, MW] = refpropm('PM','T',T_tg,'Q',0.5,'N2O');
P_sat = P_sat*1e3;
R = 8314.5/MW;

h_l_sat = refpropmJEZ('H','P',P/1e3,'Q',0,'N2O');

h_l_to_sat = h_l_sat - h_l;

% heat of vaporization (at saturation)
h_lv = h_tg_sat - h_l_sat;

% mass flow rate out via injector
% mdot_out = Cd*A_inj*sqrt(2*rho_l*(P-Po));
mdot_out = A_inj*Cd*dyer_flow(Po, P, T_l, rho_l);

% HT to saturated surface from liquid
Qdot_ls = E*Qdot('ls',T_l,T_s,rho_l,m_l,D);

% HT to gas from saturated surface
Qdot_gs = Qdot('gs',T_tg,T_s,rho_tg,m_tg,D);

% MT from surface to gas
mdot_VLB = (Qdot_ls - Qdot_gs)/(h_lv + h_l_to_sat);

% MT from gas to liquid
mdot_VLC = (P > P_sat)*(P - P_sat)*V_tg/(R*T_tg*0.1);

% if h < 1e-3
%     mdot_VLC = mdot_VLC*h/(1e-3);
% end

if mdot_VLC > 0
    %     disp(['condensation, mdot/m = ' num2str(mdot_VLC/m_tg)])
end

% enthalpy of condensed vapor entering the liquid
h_VLC = h_l_sat;

% enthalpy released when vapor condenses
dh_VLC = h_tg - h_l_sat;

% net rate of change of gas mass
mdot_tg = mdot_VLB - mdot_VLC;

% net rate of change of liquid mass
mdot_l = - mdot_VLB - mdot_out + mdot_VLC;

% HT from wall to liquid
Qdot_lw = Qdot('lw',T_l,T_lw,rho_l,m_l,D);

% net HT into liquid
Qdot_l = - Qdot_ls + Qdot_lw;

% HT into gas from wall
Qdot_gw = Qdot('gw',T_tg,T_gw,rho_tg,m_tg,D);

% net HT into gas
Qdot_tg = Qdot_gs + Qdot_gw;

% this isn't actually Udot, but Udot without the P*Vdot term (hence the i)

% there's a dh term for condensation only becouse Qdot_l has the enthalpy
% of vaporization wrapped up in it

Udot_li = - mdot_out*h_l - mdot_VLB*h_l + mdot_VLC*h_VLC + Qdot_l;

Udot_tgi = mdot_VLB*h_tg_sat - mdot_VLC*h_tg + mdot_VLC*dh_VLC + Qdot_tg;

du_drho_tg = dh_drho_tg + P/rho_tg^2  - 1/rho_tg * dP_drho_tg;

du_drho_l = dh_drho_l + P/rho_l^2  - 1/rho_l * dP_drho_l;

Vdot_l = solve_for_Vdot(Udot_tgi, u_tg, mdot_tg, m_tg, du_drho_tg, ...
    Cv_tg, Udot_li, u_l, mdot_l, m_l, du_drho_l, Cv_l, dP_drho_l, V_l, ...
    dP_dT_l, dP_dT_tg, dP_drho_tg, V_tg, P, Vdot_guess);

Vdot_tg = -Vdot_l;

Udot_tg = Udot_tgi - P*Vdot_tg;

Udot_l = Udot_li - P*Vdot_l;

rhodot_tg = mdot_tg/V_tg - m_tg/V_tg^2 * Vdot_tg;

rhodot_l = mdot_l/V_l - m_l/V_l^2 * Vdot_l;

Tdot_tg = ( ( Udot_tg - u_tg*mdot_tg )/m_tg - du_drho_tg*rhodot_tg )/Cv_tg;

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
% 2 = T_tg
% 3 = T_gw
% 4 = m_l
% 5 = T_l
% 6 = T_lw
dy = [mdot_tg; Tdot_tg; Tdot_gw; mdot_l; Tdot_l; Tdot_lw];





