function Vdot_l = solve_for_Vdot(Udot_tgi, u_tg, mdot_tg, m_tg, du_drho_tg, ...
    Cv_tg, Udot_li, u_l, mdot_l, m_l, du_drho_l, Cv_l, dP_drho_l, V_l, ...
    dP_dT_l, dP_dT_tg, dP_drho_tg, V_tg, P, Vdot_guess, Vdot_bub)
% solves for Vdot_l by using the constraint that Pdot and P for both liquid and
% vapor volumes have to be equal (eg eqns 30-37 from my 2013 JPC paper).

Vdot_l = fzero(@(Vdot_l) Vdot_eqns(Vdot_l, Udot_tgi, u_tg, mdot_tg, m_tg, du_drho_tg, ...
    Cv_tg, Udot_li, u_l, mdot_l, m_l, du_drho_l, Cv_l, dP_drho_l, V_l, ...
    dP_dT_l, dP_dT_tg, dP_drho_tg, V_tg, P, Vdot_guess, Vdot_bub), Vdot_guess, Vdot_bub);

function F = Vdot_eqns(Vdot_l, Udot_tgi, u_tg, mdot_tg, m_tg, du_drho_tg, ...
    Cv_tg, Udot_li, u_l, mdot_l, m_l, du_drho_l, Cv_l, dP_drho_l, V_l, ...
    dP_dT_l, dP_dT_tg, dP_drho_tg, V_tg, P, Vdot_guess, Vdot_bub)

Vdot_tg = -Vdot_l - Vdot_bub;

rhodot_tg = mdot_tg/V_tg - m_tg/V_tg^2 * Vdot_tg;

rhodot_l = mdot_l/V_l - m_l/V_l^2 * Vdot_l;

Udot_tg = Udot_tgi - P*Vdot_tg;

Udot_l = Udot_li - P*Vdot_l;

Tdot_tg = ( ( Udot_tg - u_tg*mdot_tg )/m_tg - du_drho_tg*rhodot_tg )/Cv_tg;

Tdot_l = ( ( Udot_l - u_l*mdot_l )/m_l - du_drho_l*rhodot_l )/Cv_l;

Vdot_l_eqn = ( dP_drho_l * (mdot_l/V_l) + dP_dT_l*Tdot_l - dP_dT_tg * Tdot_tg - dP_drho_tg * (mdot_tg/V_tg) )...
    /( dP_drho_l * (m_l/V_l^2) + dP_drho_tg * (m_tg/V_tg)^2 );

F = (Vdot_l - Vdot_l_eqn)/(abs(Vdot_guess) + 1e-6);

if isnan(F) || ~isreal(F)
    disp(['problem in vdot solution'])
end