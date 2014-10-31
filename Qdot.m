
function Q = Qdot(interface,T_fluid,T_surface,rho,m,D)
% sign convention for heat transfer: everything points towards the gas

% agw: HT from air to the gas side wall
% alw: HT from air to the liquid side wall
% ls: HT from liquid to the saturated surface
% gs: HT from the saturated surface to the gas
% lw: HT from the liquid side wall to the liquid
% gw: HT from the gas side wall to the gas

g = 9.81;

switch interface
    case {'agw','alw'}
        % from the air to the gas side wall or
        % from the air to the liquid side wall
        c = 0.59;
        n = 1/4;
        dT = T_fluid - T_surface;
        V = m/rho;
        L = V/(pi*D^2/4);
        A = pi*D*L;
        rho = 1.2;
        mu = 1.9e-5;
        k = 0.026;
        beta = 1/T_fluid;
        Cp = 1007;
        
        Ra = Cp*rho^2*g*beta*abs(dT)*L^3/(mu*k);
        
        Nu = c*Ra^n;
        
        if Ra > 1e9
            %             disp('Ra too large for AGW/ALW HT correlation')
        end
    case {'ls','gs'}
        % from the liquid to the surface or
        % from the surface to the gas
        
        T_film = .5*(T_fluid + T_surface);
        [Cp, beta, mu, k] = refpropm('CBVL','T',T_film,'D&',rho,'N2O');
        c = 0.15;
        n = 1/3;
        dT = T_fluid - T_surface;
        L = pi*D/4;
        A = pi*D^2/4;
        
        if strcmp(interface,'gs')
            dT = T_surface - T_fluid;
        end
        
        Ra = Cp*rho^2*g*beta*abs(dT)*L^3/(mu*k);
        
        Nu = c*Ra^n;
        
        if Ra > 1e11
            %             disp('Ra too large for LS/GS HT correlation')
        end
        
    case {'lw','gw'}
        % from the wall to the liquid
        % or from the wall to the gas
        T_film = 0.5*(T_fluid + T_surface);
        
        [Cp, beta, mu, k, Pr] = refpropm('CBVL^','T',T_film,'D&',rho,'N2O');
        
        c = 0.021;
        n = 2/5;
        dT = T_surface - T_fluid;
        V = m/rho;
        L = V/(pi*D^2/4);
        A = pi*D*L;
        
        Ra = Cp*rho^2*g*beta*abs(dT)*L^3/(mu*k);
        Nu = c*Ra^n;
        %         Nu = ( 0.825 + ( 0.387*Ra^(1/6) )/( 1 + ( 0.492/Pr)^(9/16) )^(8/27) )^2;
        
        
end

htc = Nu*k/L;

Q = htc*A*(dT);