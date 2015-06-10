% function G_out = injector_flow(Po, P, T_liq, rho_liq, Psat_liq, s_liq, h_liq, fluid)

function G_out = injector_flow(x_out, P, hesson_fit, rho_out, h_out, fluid, constants)

% load hesson_fit

f_feedline = constants.f_feedline;
L_feedline = constants.L_feedline;
D_feedline = constants.D_feedline;
A_inj = constants.A_inj;
Cd = constants.Cd;

G_guess = hesson_fit.fitresult(P, x_out);

G = fzero(@(G) error_fn(G), G_guess);


    function E = error_fn(G)
        u_feedline = G*Cd*A_inj / (rho_out * pi/4* D_feedline^2);

        dP_feedline = 0.5 * rho_out * u_feedline^2 * f_feedline * L_feedline/D_feedline;

        P_inj = P - dP_feedline;

        x_inj = refpropm('Q','H',h_out,'P',(P - dP_feedline)/1e3,fluid);

        G_out = hesson_fit.fitresult(P_inj, x_inj);
        
        E = (G - G_out)./G_out;

    end

end

% dyer flow model

% [rho_o, h_o] = refpropm('DH','P',6*Po/1e3,'S',s_liq,fluid);
% 
% G_HEM = rho_o*sqrt(2*(h_liq - h_o));
% 
% G_SPI = sqrt(2*rho_liq*(P - Po));
% 
% k = sqrt( (P - Po)/(Psat_liq - Po) );
% 
% G = 1/(1+k)*(k*G_SPI + G_HEM);

% % henry-fauske




% function varargout = henry_fauske_eqns(eta, P_o, x_o, v_go, v_lo)
% 
% 
% P_t = P_o*eta;
% 
% v_gt = v_go*eta^(-1/gamma);
% 
% alpha_o = x_o*v_go/( (1-x_o)*v_lo + x_o*v_go);
% alpha_t = x_o*v_gt/( (1-x_o)*v_lo + x_o*v_gt);
% 
% x_Et = ?????
% 
% N = x_Et/0.14;
% 
% if N > 1, N = 1; end
% 
% n = ( (1 - x)*c_l/c_pg + 1)/( (1 - x)*c_l/c_pg + 1/gamma);
% 
% beta = ( 1/n + (1 - v_lo/v_gt)*( (1 - x_o)*N*P_t/(x_o*(s_gE - s_lE))*ds_lE_dP ) ... 
%     - c_pg*(1/n - 1/gamma)/(s_go - s_lo) );
% 
% varargout{1} = eta - (( (1 - alpha_o)/alpha_o *(1 - eta) + gamma/(gamma - 1) )...
%     /(1/(2*beta*C^2*alpha_t^2) + gamma/(gamma-1) ) )^( gamma/(gamma-1) );
% 
% if nargout == 2
%     
% varargout{2} = sqrt( (1 - x_o)*v_lo*(P_o - P_t) + x_o*gamma/(gamma - 1)*(P_o*v_go - P_t*v_gt) )...
%     /(0.5*( (1-x_o)*v_lo + x_o*v_gt ) );
% 
% end