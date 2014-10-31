function T = get_T_from_m_U(m, U, V_tank, T_guess)

T = fzero(@(x) T_eqns(x,V_tank, m, U),[T_guess]);

function F = T_eqns(x,V_tank,m,U)
T = x(1);

[u_liq, rho_liq] = refpropm('UD','T',T,'Q',0,'N2O');
[u_vap, rho_vap] = refpropm('UD','T',T,'Q',1,'N2O');

x = (U/m - u_liq)/(u_vap - u_liq);

F = [(V_tank - m*( (1-x)/rho_liq + x/rho_vap))/V_tank];