% make a P-rho-T table for metastable calculations
% need to do for both phases

% what I really need to have is rho = rho(P,T) 
% for both liquid and vapor. 


% min point = triple
% max point = critical

function generate_PDT_table

NT = 205;
NP = 204;
% 
[Tmin] = refpropmJEZ('T','R',0,'',0,'N2O');
[Tmax] = refpropmJEZ('T','C',0,'',0,'N2O');
Tmax = Tmax - 1;
Tmin = Tmin + 1;
Pmin = refpropmJEZ('P','T',Tmin,'Q',0.5,'N2O');
Pmax = refpropmJEZ('P','T',Tmax,'Q',0.5,'N2O');
% [Dmin_liq,Dmin_vap] = refpropmJEZ('+-','T',Tmin,'Q',0.5,'N2O');
% [Dmax_liq,Dmax_vap] = refpropmJEZ('+-','T',Tmax,'Q',0.5,'N2O');

% T = linspace(Tmin,Tmax,NT);
% D_liq = linspace(Dmin_liq,Dmax_liq,ND);
% D_vap = linspace(Dmin_vap,Dmax_vap,ND);

% parfor i = 1:ND
%     
%     disp(num2str(i/ND))
%     
%     
%     for j = 1:NT
%     
%         P_liq(i,j) = refpropmJEZ('P','T',T(j),'D&',D_liq(i),'N2O');
%         P_vap(i,j) = refpropmJEZ('P','T',T(j),'D&',D_vap(i),'N2O');
%     end
% end

% for j = 1:NT
%         [Dsat_liq(j), Dsat_vap(j)] = refpropmJEZ('+-','T',T(j),'Q',0.5,'N2O');
% 
% end

% whos

% save PDT_table

% figure(1)
% contourf(T,D_liq,P_liq)
% hold on
% plot(T,Dsat_liq,'k','linewidth',3)
% colorbar
% xlabel('T')
% ylabel('D')
% title('liquid pressure')
% 
% figure(2)
% contourf(T,D_vap,P_vap)
% hold on
% plot(T,Dsat_vap,'k','linewidth',3)
% colorbar
% xlabel('T')
% ylabel('D')
% title('vapor pressure')

T = linspace(Tmin,Tmax,NT);
P = linspace(Pmin,Pmax,NP);
parfor i = 1:NP
    
    disp(num2str(i/NP))
    
%     Tsat(i) = refpropmJEZ('T','P',P(i),'Q',0.5,'N2O');
    for j = 1:NT

    [D_guess_liq, D_guess_vap] = refpropmJEZ('+-','T',T(j),'Q',0.5,'N2O');

    D_liq(i,j) = fzero(@(rho) get_D(rho,T(j),P(i)*1e3),D_guess_liq);
    D_vap(i,j) = fzero(@(rho) get_D(rho,T(j),P(i)*1e3),D_guess_vap);
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
% plot(Tsat,P,'k','linewidth',3)
% colorbar
% xlabel('P')
% ylabel('T')
% title('liquid density')
% 
% figure(2)
% % subplot(1,2,2)
% contourf(T,P,D_vap)
% hold on
% plot(Tsat,P,'k','linewidth',3)
% colorbar
% xlabel('P')
% ylabel('T')
% title('vapor density')

function F = get_D(rho, T, P)

[P_var] = refpropmJEZ('P','T',T,'D&',rho,'N2O');
P_var = P_var*1e3;

F = (P - P_var)/P;