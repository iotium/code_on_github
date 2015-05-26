% calculate system P based on m_tg, U_tg, m_l, T_l
% used for non-saturated liquid, saturated ullage
% should work for 0D or 1D
function P = get_P_from_mU_mT(m_tg, U_tg, m_l, T_l, V_tank, V_node, V_bubi, PDT, guesses)


test_val = eqns_to_solve(guesses.P);

if isnan(test_val) || (isinf(test_val) || ~isreal(test_val))
    % set P = pi so other code will catch that there's a problem
    P = pi;
    
else
    
    [P, ~, exitflag] = fzero(@(P) eqns_to_solve(P), guesses.P);
    
    if (exitflag ~= 1)
        disp('something wrong in solving for P')
    end
    
end

    function F = eqns_to_solve(P)
        
        %          rho_tg = refpropm('D', 'P', P/1e3, 'U', U_tg/m_tg, 'N2O');
        
        [rho_tg_l, rho_tg_v, u_tg_l, u_tg_v] = n2o_fits_for_getting_P(P);
        u_tg = U_tg/m_tg;
        x = (u_tg - u_tg_l)/(u_tg_v - u_tg_l);
        alpha = 1/( 1 + rho_tg_v/rho_tg_l * (1 - x)/x );
        rho_tg = alpha*rho_tg_v + (1 - alpha)*rho_tg_l;
        
%         rho_l = get_D_from_TP(T_l, P, guesses);
        
        rho_l = qinterp2(PDT.T, PDT.P, PDT.D_liq, T_l, P/1e3);
        
        V_tg = m_tg/rho_tg;
        V_l = m_l./rho_l;
        
        if length(V_bubi) == 1
            V_bub = V_bubi*V_l/(1 - V_bubi);
        else
            
            node_level = get_node_levels(V_l, V_bubi, V_node, guesses.node_level);
%             V_bub = sum(node_level.*V_bubi*V_node);
            V_bub = sum_over_nodes(V_bubi, node_level, V_node);
            
        end
        
        
        F = (V_tg + V_l + V_bub - V_tank)/V_tank;
        
%         if isnan(F) || ~isreal(F)
%             disp('problem in sol')
%         end
    end

end

