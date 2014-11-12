% make a P-rho-T table for metastable calculations
% need to do for both phases

% what I really need to have is rho = rho(P,T)
% for both liquid and vapor.


% min point = triple
% max point = critical

function generate_PDT_table

Tn_sat_max = 1.25;
Tn_sat_min = 0.8;

NT = 205;
NP = 204;
%
[Tmin, MW] = refpropm('TM','R',0,'',0,'N2O');
[Tmax] = refpropm('T','C',0,'',0,'N2O');
Tmax = Tmax - 1;
Tmin = Tmin + 1;
Pmin = refpropm('P','T',Tmin,'Q',0.5,'N2O');
Pmax = refpropm('P','T',Tmax,'Q',0.5,'N2O');

Dmax = 28.120; % [mol/L] (N2O)
Dmax = Dmax * MW; % [g/L]
Dmax = Dmax * 1; % [kg/m^3]

T = linspace(Tmin,Tmax,NT);
P = linspace(Pmin,Pmax,NP);

for i = 1:NP
    
    disp(num2str(i/NP))
    
    [D_sat_liq, D_sat_vap, Tsat(i)] = refpropm('+-T','P',P(i),'Q',0.5,'N2O');
    
    [~, ind_sat] = min(abs(T - Tsat(i)));
    
    D_liq(i,1:NT) = nan(1,NT);
    D_vap(i,1:NT) = nan(1,NT);
    
    if ind_sat > 1 && ind_sat < NT
        
        if T(ind_sat) > Tsat
            ind_sat1 = ind_sat - 1;
            ind_sat2 = ind_sat;
        else
            ind_sat1 = ind_sat;
            ind_sat2 = ind_sat + 1;
            
        end
        
        [~, ind_T_max] = min(abs(T - Tsat(i) * Tn_sat_max));
        [~, ind_T_min] = min(abs(T - Tsat(i) * Tn_sat_min));
        

        % increasing T from sat
        % expect rho to decrease
        for j = ind_sat2 : min([NT ind_T_max]) 
            if j == ind_sat2
                D_guess_liq = D_sat_liq;
                D_guess_vap = D_sat_vap;
            else
                D_guess_liq = D_liq(i,j-1);
                D_guess_vap = D_vap(i,j-1);
            end
            
            
            D_liq(i,j) = get_rho(T(j), P(i), D_guess_liq*0.97);
            D_vap(i,j) = get_rho(T(j), P(i), D_guess_vap*0.97);
            
            
        end
        
        % decreasing T from sat
        % expect rho to increase
        for j = ind_sat1:-1: max( [ 1 ind_T_min] )
            if j == ind_sat1
                D_guess_liq = D_sat_liq;
                D_guess_vap = D_sat_vap;
            else
                D_guess_liq = D_liq(i,j+1);
                D_guess_vap = D_vap(i,j+1);
            end
            
            D_liq(i,j) = get_rho(T(j), P(i), D_guess_liq/0.97);
            D_vap(i,j) = get_rho(T(j), P(i), D_guess_vap/0.97);
            
        end
        
    end
    
    %     for j = 1:NT
    %
    %         if j == 1
    %             [D_guess_liq, D_guess_vap] = refpropm('+-','T',T(j),'Q',0.5,'N2O');
    %         else
    %             D_guess_liq = D_liq(i,j-1);
    %             D_guess_vap = D_vap(i,j-1);
    %         end
    %
    % %     D_liq(i,j) = fzero(@(rho) get_D(rho,T(j),P(i)),D_guess_liq);
    % %     D_vap(i,j) = fzero(@(rho) get_D(rho,T(j),P(i)),D_guess_vap);
    %     D_liq(i,j) = fminbnd(@(rho) abs(get_D(rho,T(j),P(i))),D_guess_liq);
    %     D_vap(i,j) = fminbnd(@(rho) abs(get_D(rho,T(j),P(i))),D_guess_vap);
    %     end
end

PDT.T = T;
PDT.P = P;
PDT.D_liq = D_liq;
PDT.D_vap = D_vap;

% save('PDT_table.mat','PDT')
%
figure(1)
% subplot(1,2,1)
contourf(T,P,D_liq)
hold on
plot(Tsat,P,'k','linewidth',3)
colorbar
ylabel('P')
xlabel('T')
title('liquid density')
hold on
plot(Tsat, P, 'k-')

figure(2)
% subplot(1,2,2)
contourf(T,P,D_vap)
hold on
plot(Tsat,P,'k','linewidth',3)
colorbar
ylabel('P')
xlabel('T')
title('vapor density')
plot(Tsat, P, 'k-')


    function D = get_rho(T, P, D_guess)

%         D = fzero(@(rho) eqns_to_solve(rho), D_guess);
        D = fminbnd(@(rho) abs( eqns_to_solve(rho)), D_guess*0.8, D_guess*1.2);
        
        function F = eqns_to_solve(rho)
            
            if rho > Dmax
%                 disp('density is too high')
                delta_D = (rho - Dmax)/Dmax;
                [P_var] = refpropm('P','T',T,'D&',Dmax*0.999,'N2O');
                
                F = (P - P_var)/P;
                
                if F > 0
                    F = F + delta_D;
                else
                    F = F - delta_D;
                end
                
            else
            
                [P_var] = refpropm('P','T',T,'D&',rho,'N2O');
                F = (P - P_var)/P;
            
            end
            
            
        end
        
    end

end