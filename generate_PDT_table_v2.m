function generate_PDT_table_v2

NT = 1005;
NP = 1004;
%
[Tmin, MW] = refpropm('TM','R',0,'',0,'N2O');
[Tmax, D_c] = refpropm('TD','C',0,'',0,'N2O');
Tmax = Tmax - 1;
Tmin = Tmin + 1;
Pmin = refpropm('P','T',Tmin,'Q',0.5,'N2O');
Pmax = refpropm('P','T',Tmax,'Q',0.5,'N2O');

Dmax = 28.120; % [mol/L] (N2O)
Dmax = Dmax * MW; % [g/L]
Dmax = Dmax * 1; % [kg/m^3]

% min D is at highest T and lowest P
Dmin = refpropm('D', 'T', Tmax, 'P', Pmin, 'N2O');

T = linspace(Tmin,Tmax,NT);
P = linspace(Pmin,Pmax,NP);

ND = NP;

D_liq = nan(NP, NT);
D_vap = nan(NP, NT);

% loop through temperature
for j = 1:NT
    % find spinodal points
    
    j/NT
    
    % find the saturation P and D for a given T
    [Psat_L(j), Dsat_L(j)] = refpropm('PD','T',T(j),'Q',0,'N2O');
    [Psat_V(j), Dsat_V(j)] = refpropm('PD','T',T(j),'Q',1,'N2O');
    
    % find the D and P at the spinodal for a given T
    % for liquid it's a min of p, while for vapor it's a max
    [Dspin_L(j), Pspin_L(j)] = fminbnd(@(D) refpropm('P','T',T(j),'D&',D,'N2O'), D_c, Dsat_L(j));
    [Dspin_V(j), Pspin_V(j)] = fminbnd(@(D) -refpropm('P','T',T(j),'D&',D,'N2O'), Dsat_V(j), D_c);
    Pspin_V(j) = - Pspin_V(j); % output of fminbnd is fval, which here is -refpropm
    
    D_L = linspace(Dspin_L(j), Dmax*0.99, ND);
    D_V = linspace(Dmin, Dspin_V(j), ND);
    
    % now get P = P(rho)
    
    for k = 1:ND
        P_L(k) = refpropm('P', 'T', T(j), 'D&', D_L(k), 'N2O');
        P_V(k) = refpropm('P', 'T', T(j), 'D&', D_V(k), 'N2O');
    end
      
    % now invert it to rho = rho(P)
    
    for i = 1:NP
        % only bother if not past the spinodal
        % for liquid this means P > Pspin
        % for vapor this means P < Pspin
        
        if P(i) > Pspin_L(j)
            D_liq(i,j) = interp1(P_L, D_L, P(i),'pchip');
        end
        
        if P(i) < Pspin_V(j)
            D_vap(i,j) = interp1(P_V, D_V, P(i),'pchip');
        end
    end
    
end

PDT.T = T;
PDT.P = P;
PDT.D_liq = D_liq;
PDT.D_vap = D_vap;

save('PDT_table.mat','PDT')
%
% figure(1)
% % subplot(1,2,1)
% contourf(T,P,D_liq)
% hold on
% colorbar
% ylabel('P')
% xlabel('T')
% title('liquid density')
% hold on
% plot(T, Psat_L, 'k-')
% 
% figure(2)
% % subplot(1,2,2)
% contourf(T,P,D_vap)
% hold on
% colorbar
% ylabel('P')
% xlabel('T')
% title('vapor density')
% plot(T, Psat_L, 'k-')
        