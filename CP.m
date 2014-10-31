% casalino and pastrone model
% function [P_LRO, T_LRO, t_LRO] = CP_model(A_inj, specified_case)
function varargout = CP(varargin)

if nargin == 0
%    A_inj = 2e-5;
    specified_case = 2;
    if ~exist('A_inj','var')
        switch specified_case
            case 1
                A_inj = 2.211e-4;
            case 2
                A_inj = 2.36e-5;
            case 3
                A_inj = 7.41e-5;
            case 4
                A_inj = 7.20e-5;
            case 5
                A_inj = 4.703e-7;
            case 6
                A_inj = 2.311e-6;
        end
    end
else
    A_inj = varargin{1};
    specified_case = varargin{2};
    k_lim = varargin{3};
end

% function CP_model
% A_inj = 15e-5;
% specified_case = 1;
% problem parameters
% initialize program parameters
h = 1e-9;           % [s] initial time step
running = 1;        % [] switch, 1 = program running, 0 = program stopped
rel_tol = 1e-10;     % [] max relative error allowed in adaptive scheme
abs_tol = 1e3;     % [] max absolute error allowed in adaptive scheme
min_error = 1e-2;   % [] min error (relative to error_tol) before step size is increased
adaptive = 1;       % [] switch, 1 = use 4/5 Runge-Kutta_Fehlberg, 0 = plain 4th Runge-Kutta
h_max = 1e-0;       % [s] max allowable time step
h_min = 1e-16;      % [s] min allowable time step
t_end = 300;         % [s] end time (if LRO doesn't happen first)
LRO_tol = 1e-5;     % [s] tolerance for resolving the LRO point
non_adaptive_scheme = 0; % [] switch, 1 = 4th order Runge-Kutta, 0 = 1st order Euler
boil_tol = 1e-3;

N_dim = 6;

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

[rho_liq, rho_vap, P] = refpropmJEZ('+-P','T',Ti,'Q',0.5,'N2O');
P = P*1e3;

T_sat = Ti;

V_l = fill_level*V_tank;
V_tg = (1 - fill_level)*V_tank;

m_liq = fill_level*V_tank*rho_liq;
m_vap = (1 - fill_level)*V_tank*rho_vap;

T_liq = Ti;
T_vap = Ti;
T_lw = Ti;
T_gw = Ti;

y(:,1) = [m_vap; Ti; Ti; m_liq; Ti; Ti];

% 1 = m_tg
% 2 = T_tg
% 3 = T_gw
% 4 = m_l
% 5 = T_l
% 6 = T_lw

constants = [0; D; t_w; rho_w;
    cv_w; Cd; A_inj; Po;
    T_air; V_tank; 0; 0;
    k_w; k_lim];

% P_guess = guesses(1);
% T_tg_guess = guesses(2);
% rho_tg_guess = guesses(3);
% T_l_guess = guesses(4);
% rho_l_guess = guesses(5);

guesses = [P; rho_vap; rho_liq; 0];
P_guess = P;

derivatives = zeros(4,1);

slope = -fill_level/3;

P_lim = 0.5*(P + 1.98e5*Ti - 540e5);

P_sat_vap = P;

% begin looping
while running == 1;
    
    starti = 1;
    if n > 1
        starti = n-1;
        if n > 2
            starti = n-2;
        end
    end
    
    %     mdot_vap = bdiff(y(1,:),starti,n,t,adaptive);
    %     Tdot_vap = bdiff(y(2,:),starti,n,t,adaptive);
    %     mdot_liq = bdiff(y(4,:),starti,n,t,adaptive);
    %     Tdot_liq = bdiff(y(5,:),starti,n,t,adaptive);
    %
    if n == 1
        mdot_vap = 1e-3;
        Tdot_vap = -1;
        mdot_liq = -1;
        Tdot_liq = -1;
    else
        mdot_vap = f(1);
        Tdot_vap = f(2);
        mdot_liq = f(4);
        Tdot_liq = f(5);
    end
    
    derivatives = [mdot_vap; Tdot_vap; mdot_liq; Tdot_liq];% + 0.95*derivatives;
    
    Qdot_lw = Qdot('lw',T_liq(n),T_lw,rho_liq, m_liq(n),D);

    
    if nargin == 0
    fprintf(['t = %8.6g, P = %6.5g, T_l = %6.5g, T_tg = %6.5g, m_l = %6.5g,'...
        'm_tg = %6.5g, fill_level = %6.5g, P_lim = %6.5g, Tdot_l  = %6.5g, '...
        'm_l*Tdot_l = %6.5g, T_lw = %6.5g, T_gw = %6.5g, Qdot_lw/m_l = %6.5g\n'],...
        t(n), P(n)/6894.8, T_liq(n), T_vap(n), m_liq(n), ...
        m_vap(n), fill_level(n), P_lim(n)/6894.8, Tdot_liq, ...
        Tdot_liq*m_liq(n), T_lw, T_gw, Qdot_lw/m_liq(n));
    end
    % if we're not on the first step and error is plenty small, increase
    % step size
    if  n > 1 && adaptive == 1
        
        % if error is < min_error
        if max(abs_err/abs_tol,rel_err/rel_tol) < min_error
            
            % make h bigger
            h = min(4*h,h_max);
            
        end
        
        %         slope = (fill_level(n) - fill_level(n-1))/(t(n) - t(n-1));
%         slope = 0.05*bdiff(fill_level,starti,n,t,adaptive) + 0.95*slope
        slope_LRO = mdot_liq/V_tank/rho_liq;
        
        t_LRO = -fill_level(n)/slope_LRO + t(n);
        
        h_LRO = t_LRO - t(n); % distance to t_LRO
        
        % if the step we're about to take is >3/4 the distance to LRO
        % and the distance te LRO is bigger than the tolerance
        if (h > 0.75*h_LRO && h_LRO > LRO_tol) && (h_LRO > 0)
            
            % set h to 1/2 the distance to LRO (ie refine)
            h = 0.5*h_LRO;
            
        end
        
        P_slope = bdiff((P_sat_vap - P_lim),1,n,t,adaptive);
        t_Pswitch = -(P_sat_vap(n) - P_lim(n))/P_slope;
        h_switch = t_Pswitch - t(n);
        
        if (h > 0.75*h_switch) && (h_switch > 0 &&(h_switch > boil_tol) )
            h = h_switch/2;
        end
        
%                 T_vap_slope = bdiff(T_vap,1,n,t,adaptive);
%         
%         t_vap_low = (T_min - T_vap(n))/T_vap_slope;
%         
%         T_liq_slope = bdiff(T_liq,1,n,t,adaptive);
%         
%         t_liq_low = (T_min - T_liq(n))/T_liq_slope;
%         
%         h_T_low = min([t_vap_low - t(n); t_liq_low - t(n)]);
%         
%         
%         if (h > 0.75*h_T_low) && (h_T_low > 0)
%             h = h_T_low/2;
%         end
        
    end
    
    error_OK = 0;
    
    while error_OK == 0
        % solving differential equations
        % i = counter for
        
        constants(11) = h;
        
        for i = 1:s
            
            if i == 1
                
                % f = f( t(n) , y(n) )
                
                f = diffeqns(y(:,n), constants, derivatives, guesses, fsolve_options);
            else
                
                % f for k(2) = f( t(n) + c(2)*h , y(n) + a(2,1)*k(1) )
                % f for k(3) = f( t(n) + c(3)*h , y(n) + a(3,1)*k(1) + a(3,2)*k(2) )
                % and so on
                
                f = diffeqns(y(:,n) + sum( (ones(N_dim,1)*a(i,1:i-1)).*k(:,1:i-1),2 ), ...
                    constants, derivatives, guesses, fsolve_options);
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
                
                if (n > 1 && h_LRO < LRO_tol) && (fill_level(n) < 0.01)
                    running = 0;
                    disp('reached LRO')
                end
                
                if h < 2*h_min
                    disp('h too small, errors!')
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
    
    m_vap(n+1) = y(1,n+1);
    T_vap(n+1) = y(2,n+1);
    T_gw = y(3,n+1);
    m_liq(n+1) = y(4,n+1);
    T_liq(n+1) = y(5,n+1);
    T_lw = y(6,n+1);
    
    [rho_liq] = refpropmJEZ('+','T',T_liq(n+1),'Q',0,'N2O');
    [P_sat_liq] = refpropmJEZ('P','T',T_liq(n+1),'Q',1,'N2O');
    [rho_vap, P_sat_vap(n+1)] = refpropmJEZ('DP','T',T_vap(n+1),'Q',1,'N2O');
    
    P_sat_liq = P_sat_liq*1e3;
    P_sat_vap(n+1) = P_sat_vap(n+1)*1e3;
    
    P_lim(n+1) = 0.5*(P_sat_liq + 1.98e5*T_liq(n+1) - 540e5);
    boiling = constants(12);
    if ~boiling
        P(n+1) = P_sat_vap(n+1);
        
        if P_lim(n+1) >= P_sat_vap(n+1)
            boiling = 1;
            constants(12) = 1;
            P(n+1) = P_lim(n+1);
            disp('switched')
        end
    else
        P(n+1) = P_lim(n+1);
    end
    
    
    % if P_lim < P_sat_tg
    %     P(n+1) = P_sat_tg;
    % else
    %     P(n+1) = P_lim;
    % end
    
    fill_level(n+1) = m_liq(n+1)/rho_liq/V_tank;
    
    
    t(n+1) = t(n) + h;
    
    n = n + 1;
    
    if (sum(abs(imag(y(:,n)))) > 0) || (sum(isnan(y(:,n))) > 0)
        running = 0;
        disp('imaginary or nans')
    end
    
    if (t(n) > t_end) || ( fill_level(n) < 1e-5)
        disp('reached end t or ran out of liquid')
        running = 0;
    end
    
    %     if t(n) > ti
    %         ti = ti + 0.1;
    %         disp(num2str(t(n)))
    %     end
    
    if P(end) < 2e5
        running = 0;
        disp('pressure got low')
    end
    
    
end


if nargout > 0
    P_LRO = P(end);
    T_LRO = T_liq(end);
    t_LRO = t(end);
    varargout{1} = P_LRO/P(1);
    varargout{2} = T_LRO/T_liq(1);
    varargout{3} = t_LRO;
    
if nargout > 3
        varargout{4} = P;
        varargout{5} = T_liq;
        varargout{6} = t;
    end
else
    
    
% 1 = m_tg
% 2 = T_tg
% 3 = T_gw
% 4 = m_l
% 5 = T_l
% 6 = T_lw
    
     figure(1)
    hold on
    plot(t,P/1e6,'k-')
    xlabel('Time [s]')
    ylabel('Pressure [MPa]')
    
    figure(2)
    hold on
    plot(t,y(5,:),'k-',t,y(2,:))
    ylabel('Temperature')
    xlabel('Time [s]')
    
    legend('Liquid','Vapor')
    
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
    plot(t,T_liq,'k',t,T_vap,'b--')
    xlabel('Time [s]')
    ylabel('Temperature [K]')
    legend('Liquid','Vapor')
    
    subplot(1,3,3)
    hold on
    plot(t,m_liq,'k-',t,m_vap,'b--')
    xlabel('Time [s]')
    ylabel('Mass [kg]')
    legend('Liquid','Vapor')
% figure(1)
% hold on
% plot(t,P/6894.75729,'k-')
% title('pressure')
%
% disp('done')
end
function dy = diffeqns(y, constants, derivatives, guesses, fsolve_options)

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
boiling = constants(12);
k_w = constants(13);
k_lim = constants(14);

% retrieve derivatives calculated with backwards differencing

% retrieve variables
% 1 = m_tg
% 2 = T_tg
% 3 = T_gw
% 4 = m_l
% 5 = T_l
% 6 = T_lw
m_vap = y(1);
T_vap = y(2);
T_gw = y(3);
m_liq = y(4);
T_liq = y(5);
T_lw = y(6);

mdot_vap_guess = derivatives(1);
Tdot_vap_guess = derivatives(2);
mdot_liq_guess = derivatives(3);
Tdot_liq_guess = derivatives(4);

[rho_liq, h_liq, cp_liq, cv_liq, drho_liqdT_P] = refpropmJEZ('+HCOW','T',T_liq,'Q',0,'N2O');
[h_vap_liq, P_sat_liq] = refpropmJEZ('HP','T',T_liq,'Q',1,'N2O');
[rho_vap, h_vap, u_vap, P_sat_vap, Z_vap] = ...
    refpropmJEZ('DHUPZ','T',T_vap,'Q',1,'N2O');
h_liq_vap = refpropmJEZ('H','T',T_vap,'Q',0,'N2O');

% mult = 1.0;
% 
% c_liq = mult*c_liq;
% h_liq = mult*h_liq;
% h_liq_vap = mult*h_liq_vap;
% h_vap_liq = mult*h_vap_liq;
% h_vap = mult*h_vap;
% u_vap = mult*u_vap;

% h_liq = h_liq - 469.5e3;
% h_vap_liq = h_vap_liq - 469.5e3;
% h_vap = h_vap - 469.5e3;
% h_liq_vap = h_liq_vap - 469.5e3;
% u_vap = u_vap - 469.5e3;

P_sat_liq = P_sat_liq*1e3;
P_sat_vap = P_sat_vap*1e3;

V_vap = m_vap/rho_vap;
V_liq = m_liq/rho_liq;

P_lim = k_lim*(P_sat_liq + 1.98e5*T_liq - 540e5);

if boiling
    P = P_lim;
    
    if mdot_vap_guess < 0
        mdot_vap_guess = -10*mdot_vap_guess;
    end
else
    P = P_sat_vap;
end

% mass flow rate out via injector
mdot_out = A_inj*Cd*dyer_flow(Po, P, T_liq, rho_liq);
% mdot_out = A_inj*Cd*sqrt(2*rho_liq*(P-Po));

% HT from wall to liquid
Qdot_lw = 0;%Qdot('lw',T_liq,T_lw,rho_liq,m_liq,D);

% net HT into liquid
Qdot_liq = Qdot_lw;

% HT into gas from wall
Qdot_gw = 0;%Qdot('gw',T_vap,T_gw,rho_vap,m_vap,D);

% net HT into gas
Qdot_vap =  Qdot_gw;

delta_h_vap = h_vap - h_liq_vap;
delta_h_liq = h_vap_liq - h_liq;

dT = 0.1;
T_liq_cd = [(T_liq - 2*dT);
    (T_liq - dT);
    (T_liq + dT);
    (T_liq + 2*dT)];
T_vap_cd = T_liq_cd - T_liq + T_vap;

for j = 1:4
    [P_vap_cd(j), u_vap_cd(j), Z_cd(j)] = ...
        refpropmJEZ('PUZ','T',T_vap_cd(j),'Q',1,'N2O');
    [P_liq_cd(j), rho_liq_cd(j), u_liq_cd(j)] = ...
        refpropmJEZ('P+U','T',T_liq_cd(j),'Q',0,'N2O');
end

dP_vapdT_sat = 1e3*(1/12/dT)*(P_vap_cd(1) -8*P_vap_cd(2) + ...
    8*P_vap_cd(3) - P_vap_cd(4));
dP_liqdT_sat = 1e3*(1/12/dT)*(P_liq_cd(1) -8*P_liq_cd(2) + ...
    8*P_liq_cd(3) - P_liq_cd(4));
du_vapdT_sat = (1/12/dT)*(u_vap_cd(1) -8*u_vap_cd(2) + ...
    8*u_vap_cd(3) - u_vap_cd(4));
du_liqdT_sat = (1/12/dT)*(u_liq_cd(1) -8*u_liq_cd(2) + ...
    8*u_liq_cd(3) - u_liq_cd(4));
drho_liqdT_sat = (1/12/dT)*(rho_liq_cd(1) -8*rho_liq_cd(2) + ...
    8*rho_liq_cd(3) - rho_liq_cd(4));
dZdT_sat = (1/12/dT)*(Z_cd(1) -8*Z_cd(2) + ...
    8*Z_cd(3) - Z_cd(4));

dP_limdT = 0.5*(dP_liqdT_sat + 1.98e5);
dPdT_sat = dP_vapdT_sat;

dv_liqdT_P = -1/rho_liq^2*drho_liqdT_P;
csat_liq = cp_liq - T_liq*dv_liqdT_P*dP_liqdT_sat;

% c_liq = cp_liq;
% c_liq = cv_liq;
% c_liq = csat_liq;
c_liq = du_liqdT_sat;

IC_odes = [mdot_vap_guess, Tdot_vap_guess, mdot_liq_guess, Tdot_liq_guess];

[x, fval] = fsolve(@(x) simultaneous_system_odes(x, dZdT_sat, rho_liq, m_liq,...
    drho_liqdT_sat, dPdT_sat, dP_limdT, du_vapdT_sat, delta_h_vap, P, ...
    rho_vap, h_vap_liq, u_vap, h_liq_vap, h_liq, delta_h_liq,...
    Qdot_vap, c_liq, P_lim, Qdot_liq, mdot_out, m_vap, Z_vap, T_vap, ...
    V_vap, IC_odes, boiling), ...
    IC_odes, optimset('typicalx',IC_odes','tolfun',1e-9,'tolx',1e-9,'display','off'));

mdot_vap = x(1);
Tdot_vap = x(2);
mdot_liq = x(3);
Tdot_liq = x(4);

if fval > 1e-6
    disp('fsolve didn''t converge')
end

if  ~isreal(Tdot_liq)
    disp('problem')
end

% mass of wall exposed to liquid
m_lw = tank_wall_mass(V_liq,D,rho_w,t_w);

% mass of wall exposed to gas
m_gw = tank_wall_mass(V_vap,D,rho_w,t_w);

% HT from air to gas wall
Qdot_agw = 0;%Qdot('agw',T_air,T_gw,rho_vap,m_vap,D);

% HT from air to liquid wall
Qdot_alw = 0;%Qdot('alw',T_air,T_lw,rho_liq,m_liq,D);

Vdot_vap = -mdot_liq/rho_liq + m_liq/rho_liq^2*drho_liqdT_sat*Tdot_liq;

% rate of change of mass of gass wall
mdot_gw = 4*Vdot_vap*t_w*rho_w/D;

% conduction from liquid wall to gas wall
L_tank = 4*V_tank/(pi*D^2);
L_wc = L_tank/2;
Qdot_wc = k_w*(T_lw - T_gw)*pi*D*t_w/L_wc;

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
dy = [mdot_vap; Tdot_vap; Tdot_gw; mdot_liq; Tdot_liq; Tdot_lw];

function F = simultaneous_system_odes(x, dZdT_sat, rho_liq, m_liq,...
    drho_liqdT_sat, dPdT_sat, dP_limdT, du_vapdT_sat, delta_h_vap, P, ...
    rho_vap, h_vap_liq, u_vap, h_liq_vap, h_liq, delta_h_liq,...
    Qdot_vap, c_liq, P_lim, Qdot_liq, mdot_out, m_vap, Z_vap, T_vap, ...
    V_vap, IC, boiling)

mdot_vap = x(1);
Tdot_vap = x(2);
mdot_liq = x(3);
Tdot_liq = x(4);

Zdot_vap = dZdT_sat*Tdot_vap;

Vdot_vap = -mdot_liq/rho_liq + m_liq/rho_liq^2*drho_liqdT_sat*Tdot_liq;

if boiling
    % boiling has started
    Pdot = dP_limdT*Tdot_liq;
    mdot_cond = 0;
    mdot_evap = mdot_vap;
    
else
    % boiling hasn't started yet
    Pdot = dPdT_sat*Tdot_vap;
    mdot_cond = -mdot_vap;
    mdot_evap = 0;
    
end

% RHS_real_gas = mdot_vap/m_vap + Zdot_vap/Z_vap + Tdot_vap/T_vap;
% LHS_real_gas = Pdot/P + Vdot_vap/V_vap;
% RHS_vap_energy = -P*Vdot_vap + Qdot_vap;
% LHS_vap_energy = m_vap*du_vapdT_sat*Tdot_vap - mdot_cond*(delta_h_vap - P/rho_vap) - mdot_evap*(h_vap_liq - u_vap);
% RHS_liq_energy = mdot_cond*(h_liq_vap - h_liq) - mdot_evap*(delta_h_liq) + Qdot_liq;
% LHS_liq_energy = c_liq*m_liq*Tdot_liq;
% RHS_liq_m = mdot_liq;
% LHS_liq_m = mdot_cond - mdot_evap - mdot_out;

% F = [ (RHS_real_gas - LHS_real_gas)/(RHS_real_gas + 1e-9);
%     (RHS_vap_energy - LHS_vap_energy)/(RHS_vap_energy + 1e-9);
%     (RHS_liq_energy - LHS_liq_energy)/(RHS_liq_energy + 1e-9);
%     (RHS_liq_m - LHS_liq_m)/(RHS_liq_m + 1e-9)];

mdot_vap_real_gas = m_vap*(Pdot/P + Vdot_vap/V_vap - Zdot_vap/Z_vap - Tdot_vap/T_vap);
Tdot_vap_energy =  (-P*Vdot_vap + Qdot_vap + ...
    mdot_cond*(delta_h_vap - P/rho_vap) + mdot_evap*(h_vap_liq - u_vap) )...
    /(m_vap*du_vapdT_sat);
Tdot_liq_energy = (mdot_cond*(h_liq_vap - h_liq) - mdot_evap*(delta_h_liq) ...
    + Qdot_liq)/(c_liq*m_liq);
mdot_liq_m = mdot_cond - mdot_evap - mdot_out;

F = [ (mdot_vap_real_gas - mdot_vap);
    (Tdot_vap_energy - Tdot_vap);
    (Tdot_liq_energy - Tdot_liq);
    (mdot_liq_m - mdot_liq)];

F = F./(IC(:) + 1e-6);
