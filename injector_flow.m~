function G = injector_flow(Po, P, T_liq, rho_liq, Psat_liq, s_liq, h_liq)


% [Psat_liq, s_liq, h_liq] = refpropm('PSH','T',T_liq,'Q',0,'N2O');
% Psat_liq = 1e3*Psat_liq;

[rho_o, h_o] = refpropm('DH','P',Po/1e3,'S',s_liq,'N2O');

G_HEM = rho_o*sqrt(2*(h_liq - h_o));

G_SPI = sqrt(2*rho_liq*(P - Po));

k = sqrt( (P - Po)/(Psat_liq - Po) );

G = 1/(1+k)*(k*G_SPI + G_HEM);

