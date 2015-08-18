function Pi = converge_on_IC(Ti, V_tg, V_l, V_tank, V_node, V_bubi, ...
    PDT, guesses, constants, fluid)

% given initial temp, need to find P, rho_l, rho_tg, m_tg, U_tg, m_l, T_l
% don't want to just use saturation properties because it comes out wrong

P_guess = 1e3*refpropm('P','T',Ti,'Q',0.5,fluid);

guesses.P = P_guess;

Pi = fzero(@equations_to_solve, P_guess, constants.fsolve_options);

    function E = equations_to_solve(P1)
        
        if strcmp(constants.property_source,'PDT')
              
            rho_l = qinterp2(PDT.T, PDT.P, PDT.D_liq, Ti, P1/1e3);
            
            % if rho_l is NaN, it means we went outside the bounds of PDT, so
            % instead extrapolate it using interp2 (slower than qinterp2)
            if isnan(rho_l)
                rho_l = interp2(PDT.T, PDT.P, PDT.D_liq, Ti, P1/1e3, 'spline');
            end
            
            [~, rho_tg, ~, u_tg] = fits_for_getting_P(P1, fluid);

            
        elseif strcmp(constants.property_source,'refprop')
            
            rho_l = get_D_from_TP(Ti, P1, guesses, constants, fluid);
            
            [rho_tg, u_tg] = refpropm('DU','P',P1/1e3,'Q',1,fluid);
            
        end
                        
        
        
        rho_tg_sat = rho_tg;
        
        m_tg = rho_tg * V_tg;
        m_l = rho_l * V_l;
        
        U_tg = m_tg * u_tg;
        T_l = Ti;
        
        P2 = get_P_from_mU_mT(m_tg, U_tg, m_l, T_l, ...
            V_tank, V_node, V_bubi, fluid, PDT, constants, guesses);
        
        E = (P1 - P2)/P1;
        
    end

end