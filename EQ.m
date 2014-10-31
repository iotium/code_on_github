% function EQ_model_wall_HT_v4
% function [P_LRO, T_LRO, t_LRO] = EQ_model_wall_HT_v4(A_inj, specified_case)
function varargout = EQ(varargin)

if nargin == 0
    %     A_inj = 1.764e-4;
    specified_case = 2;
    if ~exist('A_inj','var')
        switch specified_case
            case 1
                A_inj = 1.764e-4;
            case 2
                A_inj = 1.97e-5;
            case 3
                A_inj = 6.01e-5;
            case 4
                A_inj = 6.41e-5;
            case 5
                A_inj = 3.167e-7;
            case 6
                A_inj = 1.413e-6;
        end
    end
else
    A_inj = varargin{1};
    specified_case = varargin{2};
end

% disp(num2str(A_inj))
% specified_case = 5;
% A_inj = 3.0123e-7;

% tic

% specified_case = 4; % [] 0 = no specified case, use following
% 1 = greg's setup from his JPC paper on modeling
% 2 = my tank testing experimental setup

% [m] tank diameter
% [m] tank diameter

% initialize program parameters
h = 1e-4;           % [s] initial time step
running = 1;        % [] switch, 1 = program running, 0 = program stopped
rel_tol = 1e-6;     % [] max relative error allowed in adaptive scheme
abs_tol = 1e-1;     % [] max absolute error allowed in adaptive scheme
min_error = 1e-3;   % [] min error (relative to error_tol) before step size is increased
adaptive = 1;       % [] switch, 1 = use 4/5 Runge-Kutta_Fehlberg, 0 = plain 4th Runge-Kutta
h_max = 1e-0;       % [s] max allowable time step
h_min = 1e-16;      % [s] min allowable time step
t_end = 300;         % [s] end time (if LRO doesn't happen first)
LRO_tol = 1e-6;     % [s] tolerance for resolving the LRO point
non_adaptive_scheme = 1; % [] switch, 1 = 4th order Runge-Kutta, 0 = 1st order Euler

N_dim = 4;

fsolve_options = optimset('display','off');

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
        fill_level = 0.8;        % [] initial fill_level ratio (by volume)
        %         E = 2e2;          % [] heat transfer multiplier
        V_tank = 0.05;   % [m^3] tank volume
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
    
    bs = [16/135; 0;  6656/12825;    28561/56430;    -9/50;  2/55]; % 5th order
    
    b = [25/216; 0; 1408/2565;     2197/4104;      -1/5;   0]; % 4th order
    
end

s = length(c);

k = zeros(N_dim,s);
ti = 0;
n = 1;              % [] counter

t = 0;

[rho_liq, rho_vap, P] = refpropm('+-P','T',Ti,'Q',0.5,'N2O');
P = P*1e3;

m_l = fill_level*V_tank*rho_liq;
m_tg = (1 - fill_level)*V_tank*rho_vap;

m = fill_level*V_tank*rho_liq + (1 - fill_level)*V_tank*rho_vap;

X = (1 - fill_level)*V_tank*rho_vap/m;

V_tg = (1 - fill_level)*V_tank;

[u] = refpropm('U','T',Ti,'Q',X,'N2O');

T = Ti;

U = u*m;

y(1,1) = m;
y(2,1) = U;
y(3,1) = T;
y(4,1) = T;

constants = [0; D; t_w; rho_w;
    cv_w; Cd; A_inj; Po;
    T_air; V_tank; 0; k_w];

derivatives = 0;

guesses = [Ti; fill_level];

while running == 1
    
    Vdot_tg(n+1) = V_tank*bdiff(V_tg/V_tank,1,n,t,adaptive);
    
    derivatives = Vdot_tg(n+1);
    if nargin == 0
        fprintf(['t = %6.4g, P = %6.4g, T = %6.4g, m = %6.4g,'...
            'fill_level = %6.4g, Vdot_tg = %6.4g,'...
            'u_l = %6.4g\n'],...
            t(n), P(n)/6894.8, T(n), y(1,end), fill_level(n), ...
            Vdot_tg(n+1), y(2,end)/y(1,end));
        
    end
    % if we're not on the first step and error is plenty small, increase
    % step size
    if  n > 1 && adaptive == 1
        
        % if error is < min_error
        if max(abs_err/abs_tol,rel_err/rel_tol) < min_error
            
            % make h bigger
            h = min(4*h,h_max);
            
        end
        
        slope = (fill_level(n) - fill_level(n-1))/(t(n) - t(n-1));
        
        t_LRO = -fill_level(n)/slope + t(n);
        
        h_LRO = t_LRO - t(n); % distance to t_LRO
        
        % if the step we're about to take is >3/4 the distance to LRO
        % and the distance te LRO is bigger than the tolerance
        if h > 0.75*h_LRO && h_LRO > LRO_tol
            
            % set h to 1/2 the distance to LRO (ie refine)
            h = 0.5*h_LRO;
            
        end
        
    end
    
    
    error_OK = 0;
    
    while error_OK == 0
        % solving differential equations
        % i = counter for
        
        for i = 1:s
            
            if i == 1
                
                % f = f( t(n) , y(n) )
                
                f = diffeqns(y(:,n), constants, derivatives, guesses);
            else
                
                % f for k(2) = f( t(n) + c(2)*h , y(n) + a(2,1)*k(1) )
                % f for k(3) = f( t(n) + c(3)*h , y(n) + a(3,1)*k(1) + a(3,2)*k(2) )
                % and so on
                
                f = diffeqns(y(:,n) + sum( (ones(N_dim,1)*a(i,1:i-1)).*k(:,1:i-1),2 ), ...
                    constants, derivatives, guesses);
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
        
        if adaptive == 1
            % using adaptive scheme, need to check error
            
            err = k*(b - bs);   % absolute error (diff. between 5th and 4th order estimates)
            
            rel_err = abs(err./( y(:,n) ));  % relative error
            
            rel_err = max(rel_err(isfinite(rel_err)));  % fix rel_err to the maximum finite value of rel_err
            
            abs_err = max(abs(err(isfinite(err)))); % do the same for abs_err
            
            % check for possible problems: isempty statements are in case
            % abs and rel err are both full of non-finite values
            % isnan checks for nan's
            % isreal checks for imaginary numbers
            error_conditions = isempty(rel_err) + isempty(abs_err) +  sum(isnan(rel_err + abs_err)) + sum(1 - isreal(err));
            
            % if any of those fail, set rel_err large so that the step gets
            % recomuputed
            if error_conditions > 0
                rel_err = 1;
            end
            
            
            
            if ( rel_err < rel_tol && abs_err < abs_tol) || (h < 1.25*h_min)
                % meeting the error requirement or step size is too
                % small already
                error_OK = 1;
                
                if n > 1 && h < LRO_tol
                    running = 0;
                end
                
                if h < 2*h_min
                    disp('h too small, errors!')
                end
                
            else
                
                % not meeting error requirements or LRO is close
                
                if rel_err == 0 || abs_err == 0
                    
                    sh = 0.1;
                    
                else
                    
                    if rel_err/rel_tol > abs_err/abs_tol
                        
                        sh = 0.84*( rel_tol*h / (2*rel_err) )^(1/4);
                    else
                        sh = 0.84*( abs_tol*h / (2*abs_err) )^(1/4);
                    end
                end
                
                if sh < 0.1
                    sh = 0.1;
                elseif sh > 4.0
                    sh = 4.0;
                end
                
                h = h*sh;
                
                h_min = 16*eps(t(n));
                
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
    
    T(n+1) = get_T_from_m_U(y(1,n+1), y(2,n+1), V_tank, T(n));
    
    U = y(2,n+1);
    m = y(1,n+1);
    
    [u_l, rho_l, P(n+1), h_l] = refpropm('UDPH','T',T(n+1),'Q',0,'N2O');
    [u_tg, rho_tg] = refpropm('UD','T',T(n+1),'Q',1,'N2O');
    P(n+1) = P(n+1)*1e3;
    x = (U/m - u_l)/(u_tg - u_l);
    
    m_l(n+1) = (1-x)*m;
    m_tg(n+1) = x*m;
    
    V_tg(n+1) = m_tg(n+1)/rho_tg;
    V_l(n+1) = m_l(n+1)/rho_l;
    
    fill_level(n+1) = V_l(n+1)/V_tank;
    
    guesses = [T(n+1), fill_level(n+1)];
    
    t(n+1) = t(n) + h;
    
    n = n + 1;
    
    if (sum(imag(y(:,n))) > 0) || (sum(isnan(y(:,n))) > 0)
        running = 0;
        disp('imaginary or nans')
    end
    
    
    if (t(n) > t_end) || (fill_level(n) < 1e-5)
        disp('reached end')
        running = 0;
    end
    
    if t(n) > ti
        ti = ti + 1;
        %         disp(num2str(t(n)))
    end
    
    
end

% toc

% figure(1)
% hold on
% plot(t,P/6895,'k--')
%
% figure(2)
% hold on
% plot(t,T,'r-.')
if nargout > 0
    P_LRO = P(end);
    T_LRO = T(end);
    t_LRO = t(end);
    varargout{1} = P_LRO/P(1);
    varargout{2} = T_LRO/T(1);
    varargout{3} = t_LRO;
    
    if nargout > 3
        varargout{4} = P;
        varargout{5} = T;
        varargout{6} = t;
    end
else
    figure(1)
    hold on
    plot(t,P/1e6,'k-')
    xlabel('Time [s]')
    ylabel('Pressure [MPa]')
    
    figure(2)
    hold on
    plot(t,T,'r-')
    ylabel('Temperature')
    xlabel('Time [s]')
    
    figure(3)
    hold on
    plot(t,y(4,:),'k-',t,y(3,:),'b--')
    title('wall temp')
    xlabel('Time [s]')
    
    legend('liquid','vapor')
    
    figure(4)
    hold on
    plot(t,m_l,'k-',t,m_tg,'b--')
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
    plot(t,T)
    xlabel('Time [s]')
    ylabel('Temperature [K]')
    
    subplot(1,3,3)
    hold on
    plot(t,m_l,'k-',t,m_tg,'b--')
    xlabel('Time [s]')
    ylabel('Mass [kg]')
    legend('Liquid','Vapor')
    
end

function dy = diffeqns(y, constants, derivatives, guesses)
m = y(1);
U = y(2);
T_gw = y(3);
T_lw = y(4);

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
k_w = constants(12);

Vdot_tg = derivatives;

T_guess = guesses(1);

T = get_T_from_m_U(m, U, V_tank, T_guess);

[u_l, rho_l, P, h_l] = refpropm('UDPH','T',T,'Q',0,'N2O');
[u_tg, rho_tg] = refpropm('UD','T',T,'Q',1,'N2O');

x = (U/m - u_l)/(u_tg - u_l);

P = P*1e3;

% mdot = Cd*A_inj*sqrt(2*rho_l*(P - Po));

mdot = Cd*A_inj*dyer_flow(Po, P, T, rho_l);

% need: rho_l, m_l, T_l, rho_tg, m_tg, T_tg, V_l, V_tg, Vdot_tg
T_l = T;
T_tg = T;

% rho_tg = refpropm('D','T',T,'Q',1,'N2O');
m_l = m*(1-x);
m_tg = m*x;

V_l = m_l/rho_l;
V_tg = m_tg/rho_tg;

% HT from wall to liquid
Qdot_lw = Qdot('lw',T_l,T_lw,rho_l,m_l,D);

% HT into gas from wall
Qdot_gw = Qdot('gw',T_tg,T_gw,rho_tg,m_tg,D);

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
Tdot_gw = (Qdot_agw - Qdot_gw + Qdot_wc +cv_w*mdot_gw*(T_lw - T_gw))/(m_gw*cv_w);

% rate of change of temperature of liquid wall
Tdot_lw = (Qdot_alw - Qdot_lw - Qdot_wc)/(m_lw*cv_w);

Udot = -h_l*mdot + Qdot_lw + Qdot_gw;

dy = [-mdot; Udot; Tdot_gw; Tdot_lw];

