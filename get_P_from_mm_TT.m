function P = get_P_from_mm_TT(m_tg, T_tg, m_l, T_l, V_tank, V_bubi, PDT, guesses)
% calculates pressure based on liquid and vapor masses + temperatures
% uses the constraint that the volumes have to add up to the tank volume

P_guess = guesses(1);
rho_guesses = guesses(2:3);

[P, ~, exitflag] = fzero(@(P) P_eqns(P, m_tg, T_tg, m_l, T_l, V_tank, ...
    V_bubi, PDT, rho_guesses), P_guess);

if (exitflag ~= 1)
    disp('something wrong in solving for P')
end

function F = P_eqns(P, m_tg, T_tg, m_l, T_l, V_tank, V_bubi, PDT, rho_guesses)

% rho_tg_guess = rho_guesses(1);
% rho_l_guess = rho_guesses(2);

% rho_l = easy_D(P, T_l, PDT, rho_l_guess, fsolve_options);
% rho_tg = easy_D(P, T_tg, PDT, rho_tg_guess, fsolve_options);

rho_l = interp2(PDT.T,PDT.P,PDT.D_liq,T_l,P/1e3,'linear');
rho_tg = interp2(PDT.T,PDT.P,PDT.D_vap,T_tg,P/1e3,'linear');

V_l = m_l/rho_l;

V_tg = m_tg/rho_tg;

F = (V_tg + V_l - V_tank - V_l*V_bubi)/V_tank;

if isnan(F) || ~isreal(F)
    disp('problem')
end