function Vdot_l = solve_for_Vdot(Udot_tgi, mdot_tg, m_tg, ...
    Udot_li, u_l, mdot_l, m_l, du_drho_l, Cv_l, dP_drho_l, V_l, ...
    dP_dT_l, V_tg, P, drho_dx_P_tg, drho_dP_x_tg, u_tg_v_sat, u_tg_l_sat, x, ...
    du_dT_sat_tg_v, du_dT_sat_tg_l, dP_dT_tg_sat, rho_tg_v, u_tg, guesses, Vdot_bub)
% solves for Vdot_l by using the constraint that Pdot and P for both liquid and
% vapor volumes have to be equal (eg eqns 30-37 from my 2013 JPC paper).
% this one assumes that the ullage is a saturated two phase mixture
% see page 42 of notebook #5

C4 = ( x * du_dT_sat_tg_v + (1-x)*du_dT_sat_tg_l -(u_tg_v_sat - u_tg_l_sat)/drho_dx_P_tg * drho_dP_x_tg * dP_dT_tg_sat ); 

C5 = (u_tg_v_sat - u_tg_l_sat)/drho_dx_P_tg;

Vdot_l = (V_l* ((Cv_l *dP_dT_tg_sat* m_l* (Udot_tgi - mdot_tg* u_tg + P *Vdot_bub)* V_l + ...
        C4 * m_tg* (-Cv_l* mdot_l *dP_drho_l* m_l + mdot_l *dP_dT_l *du_drho_l *m_l - ...
           dP_dT_l *Udot_li* V_l + mdot_l* dP_dT_l *u_l *V_l))* V_tg^2 - ...
     C5 *Cv_l* dP_dT_tg_sat* m_l *m_tg* V_l*(m_tg *Vdot_bub + ...
        mdot_tg* V_tg)))/(C5* Cv_l* dP_dT_tg_sat* m_l* m_tg^2 *V_l^2 - (C4* (Cv_l *...
dP_drho_l - dP_dT_l* du_drho_l) *m_l^2 *m_tg + (Cv_l* dP_dT_tg_sat *m_l + ...
         C4 *dP_dT_l* m_tg)* P*V_l^2) *V_tg^2); 

%{
% the "guess" variables use the guesses.Vdot_l value
% to get an estimate for Pdot and use that to normalize the error.
Udot_l_guess = Udot_li - P*guesses.Vdot_l;
rhodot_l_guess = mdot_l/V_l - m_l/V_l^2 * guesses.Vdot_l;
Tdot_l_guess = ( ( Udot_l_guess - u_l*mdot_l )...
    /m_l - du_drho_l*rhodot_l_guess )/Cv_l;
Pdot_l_guess = dP_dT_l*Tdot_l_guess + dP_drho_l*rhodot_l_guess;



test_val = Vdot_eqns(guesses.Vdot_l);

if isnan(test_val) || (isinf(test_val) || ~isreal(test_val))
    
    Vdot_l = pi;
    
else
    
    Vdot_l = fzero(@(Vdot_l) Vdot_eqns(Vdot_l), guesses.Vdot_l, Vdot_bub);
    
end

    function F = Vdot_eqns(Vdot_l)
        % equation that fzero tries to solve
        
        rhodot_l = mdot_l/V_l - m_l/V_l^2 * Vdot_l;
        
        Udot_l = Udot_li - P*Vdot_l;
        
        udot_l = ( Udot_l - u_l*mdot_l)/m_l;
        
        Tdot_l = ( udot_l - du_drho_l*rhodot_l )/Cv_l;
        
        Pdot_l = dP_dT_l*Tdot_l + dP_drho_l*rhodot_l;
        
        Vdot_tg = -Vdot_l - Vdot_bub;
        
        rhodot_tg = mdot_tg/V_tg - m_tg/V_tg^2 * Vdot_tg;
        
        Udot_tg = Udot_tgi - P*Vdot_tg;
        
        udot_tg = (Udot_tg - u_tg*mdot_tg)/m_tg;
        
        xdot_tg = V_tg/m_tg * rhodot_tg + rho_tg_v/m_tg * Vdot_tg - rho_tg_v*V_tg/m_tg^2 * mdot_tg ;
        
%         Tdot_tg = ( Udot_tg - (u_tg_v_sat - u_tg_l_sat) / ( V_tg * drho_dx_P_tg )...
%             * ( mdot_tg - m_tg/V_tg * Vdot_tg ) )...
%             /( x*du_dT_sat_tg_v + (1-x)*du_dT_sat_tg_l - ...
%             (u_tg_v_sat - u_tg_l_sat)/drho_dx_P_tg * drho_dP_x_tg * dP_dT_tg_sat );
        
%         Tdot_tg = ( udot_tg - (u_tg_v_sat - u_tg_l_sat)*xdot_tg )...
%             /(du_dT_sat_tg_v * x + (1 - x)*du_dT_sat_tg_l);

        Tdot_tg = (-(u_tg_v_sat - u_tg_l_sat)/drho_dx_P_tg * (mdot_tg - m_tg/V_tg*Vdot_tg)*1/V_tg + udot_tg )...
            ./( x * du_dT_sat_tg_v + (1-x)*du_dT_sat_tg_l -(u_tg_v_sat - u_tg_l_sat)/drho_dx_P_tg * drho_dP_x_tg * dP_dT_tg_sat ); 


        Pdot_tg = dP_dT_tg_sat * Tdot_tg;
        
        F = (Pdot_l - Pdot_tg)/(abs(Pdot_l_guess) + 1e-6);
        
        %         if isnan(F) || ~isreal(F)
        %             disp(['problem in vdot solution'])
        %         end
        
    end



Vdot_l_math = (V_l* ((Cv_l *dP_dT_tg_sat* m_l* (Udot_tgi - mdot_tg* u_tg + P *Vdot_bub)* V_l + ...
        C4 *m_tg* (-Cv_l* mdot_l *dP_drho_l* m_l + mdot_l *dP_dT_l *du_drho_l *m_l - ...
           dP_dT_l *Udot_li* V_l + mdot_l* dP_dT_l *u_l *V_l))* V_tg^2 - ...
     C5 *Cv_l* dP_dT_tg_sat* m_l *m_tg* V_l*(m_tg *Vdot_bub + ...
        mdot_tg* V_tg)))/(C5* Cv_l* dP_dT_tg_sat* m_l* m_tg^2 *V_l^2 - (C4* (Cv_l *...
dP_drho_l - dP_dT_l* du_drho_l) *m_l^2 *m_tg + (Cv_l* dP_dT_tg_sat *m_l + ...
         C4 *dP_dT_l* m_tg)* P*V_l^2) *V_tg^2); 
      
fprintf('comparison with mathematica: %0.6g\n',(Vdot_l - Vdot_l_math)/Vdot_l)
%}
end