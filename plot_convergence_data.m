% plotting data for thesis
% load the model data before running
close all
clear all 

% waiting on N_nodes = 256, 512;

% rel_tol
% 10^-1	18
% -2	16
% -3	19
% -5	30
% -6    28
% 
% N_mom
% 4     21
% 6     24
% 8     25
% 10    27
% 
% N_nodes
% 8     17
% 16	15
% 32	21
% 64	23
% 128   31

% N_rw
% 8     48
% 16    49
% 32    50
% 64    51  
% 128   52


% % wall points
% data_files = [48 49 50 51 52];
% N_rw = [8 16 32 64 128];
% ref_file = 52;
% convergence_x = N_rw;
% x_name = 'Number of Nodes in Wall';

% % moments
% data_files = [21 24 25 27];
% N_mom = [4 6 8 10];
% ref_file = 27;
% convergence_x = N_mom;
% x_name = 'Number of Moments';

% % nodes
% data_files = [17 15 21 23 31];% 29];
% N_nodes = [8 16 32 64 128];% 256];
% ref_file = 31;
% convergence_x = N_nodes;
% x_name = ['Number of Nodes'];

% rel tol
data_files = [18 16 19 22 30 28];
rel_tol = 10.^[-1 -2 -3 -4 -5 -6];
ref_file = 28;
convergence_x = rel_tol;
x_name = ['Relative Tolerance'];

cd('/Users/jez/School/stanford/compiled research/tank modeling/data from server v6/convergence study data files');

data = load(['bubble_sim_data' num2str(ref_file)]);

P_ref = data.P;
t_ref = data.t;

P_LL_25_ref = interp1(data.fill_level, data.P, 0.25);

P_LL_25 = nan*ones(size(data_files));

plot_str = {'k', 'b', 'r', 'k--', 'b--', 'r--'};

for i = 1:length(data_files)
    clear data
    filename = ['bubble_sim_data' num2str(data_files(i))];
    data = load(filename);
    
    P_compare = interp1(data.t, data.P, t_ref);
    
    figure(1)
    hold on
    plot(data.t, data.P/1e6, plot_str{i})
    
    figure(2)
    hold on
    plot(t_ref, abs((P_compare - P_ref))./P_ref, plot_str{i})
    set(gca,'yscale','log')
    
    ind_OK = find(~isnan(P_compare));
    
    legend_str{i} = sprintf('%0.5g',convergence_x(i));
    
    error_int(i) = trapz(t_ref(ind_OK), abs((P_compare(ind_OK) - P_ref(ind_OK))));
        
% t_LL_25(i) = interp1(data.fill_level, data.t, 0.25);

    if min(data.fill_level) < 0.05

        [~, ind_guess] = min(abs(data.fill_level - 0.25));

        P_LL_25(i) = interp1(data.fill_level(ind_guess+[-20:20]), data.P(ind_guess+[-20:20]), 0.25);

    end
end

int_ref = trapz(t_ref, P_ref);

figure(1)
xlabel('Time [s]')
ylabel('Pressure [MPa]')
legend(legend_str)
make_text_big('shift_start')
make_text_big

figure(2)
xlabel('Time [s]')
ylabel('Relative Error in P')
legend(legend_str{1:end-1})

make_text_big('shift_start')
make_text_big

% figure(3)
% plot(convergence_x, abs(P_LL_25 - P_LL_25_ref)./P_LL_25_ref, 'ks-')
% set(gca,'yscale','log')
% xlabel(x_name)
% ylabel('Relative Error in P at 25% F.L.')

% if strcmp(x_name,'Relative Tolerance')
% set(gca,'xscale','log')
% end
% make_text_big


figure(4)
plot(convergence_x(1:end-1), error_int(1:end-1)/int_ref, 'ks-')
set(gca,'yscale','log')
xlabel(x_name)
ylabel('Integrated Relative Error')


if strcmp(x_name,'Relative Tolerance')
set(gca,'xscale','log')
end
make_text_big



% 
% [P_exp, T_lw_out_exp, LL_exp] = load_experimental_data(t, specified_case);
% 
% m_out = cumtrapz(t, mdot_out_liq + mdot_out_vap);
% mh_out = cumtrapz(t, (mdot_out_liq.*h_l + mdot_out_vap.*h_tg_sat));
% 
% T_l = y(4+N_rw,:);
% m_l = y(3+N_rw,:);
% 
% m_tg = y(1,:);
% 
% T_lw_in = y(5+N_rw,:);
% T_gw_in = y(3,:);
% T_lw_out = y(4+2*N_rw,:);
% 
% % P vs t
% figure(1)
% plot(t,P/1e6,'k', t, P_exp/1e6,'b')
% legend('Model','Experiment')
% ylabel('Pressure [MPa]')
% xlabel('Time [s]')
% make_text_big('shift_start')
% 
% 
% % fill level vs t
% figure(2)
% plot(t, fill_level*100, 'k', t, LL_exp*100, 'b')
% xlabel('Time [s]')
% ylabel('Fill Level [%]')
% legend('Model','Experiment')
% make_text_big
% 
% % % find P at 25% fill level
% % t_LL_25 = interp1(fill_level, t, 0.25);
% % P_LL_25 = interp1(fill_level, P, 0.25);
% 
% figure(3)
% plot(t,T_l-273.15,'k',t,T_s-273.15,'k--')
% legend('Liquid','T_{sat}(P) = T_{tg}')
% ylabel('Temperature [C]')
% xlabel('Time [s]')
% make_text_big
% 
% figure(4)
% plot(t, mdot_out_liq + mdot_out_vap, 'k')
% xlabel('Time [s]')
% ylabel('Mass Flow Rate [kg/s]')
% make_text_big
% 
% figure(5)
% plot(t, x_tg,'k')
% xlabel('Time [s]')
% ylabel('Ullage Vapor Mass Fraction')
% make_text_big('shift_start')
% 
% % find sauter mean diameter, void fraction, and number density
% % versus space, at time indices of t_plot
% SMD = NaN*ones(N_nodes,length(t_plot));
% alpha = SMD;
% ND = SMD;
% % dx_nodes = 1/(N_nodes-1);
% x_nodes = linspace(0,1,N_nodes)';
% 
% t_plot = linspace(0,1,10)*t(end);
% 
% for i = 1:length(t_plot)
%     [~,i_plot] = min(abs(t - t_plot(i)));
%     ind_full = 1:(N_full(i_plot));
%     SMD(ind_full,i) = mom(ind_full,4,i_plot)./mom(ind_full,3,i_plot);
%     alpha(ind_full,i) = V_bubi(ind_full,i_plot);
%     ND(ind_full,i) = mom(ind_full,1,i_plot);
%     
%     figure(6)
%     hold on
%     plot(x_nodes, SMD(:,i), 'k')
% 
%     figure(7)
%     hold on
%     plot(x_nodes, alpha(:,i), 'k')
% 
%     figure(8)
%     hold on
%     plot(x_nodes, ND(:,i), 'k')
% 
% end
% 
% figure(6)
% xlabel('Normalized Height')
% ylabel('Sauter Mean Diameter [m]')
% make_text_big
% 
% figure(7)
% xlabel('Normalized Height')
% ylabel('Void Fraction')
% make_text_big
% 
% figure(8)
% xlabel('Normalized Height')
% ylabel('Bubble Number Density [#/m^3]')
% make_text_big
% 
% % 
% % figure(6)
% % plot(t, gas_holdup_injector, 'k')
% % xlabel('Time [s]')
% % ylabel('Outlet Vapor Mass Fraction')
% % make_text_big
% 
% 
% % figure(4)
% % plot(t,T_lw_in,'k-',t,T_gw_in,'b--',t,T_lw_out,'r-.',t,T_lw_out_exp+273.15,'k:')
% % title('wall temp')
% % xlabel('Time [s]')
% % legend('liquid','vapor','outside liquid','experimental outside liquid')
% % title('wall temp')