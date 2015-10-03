% plotting data for thesis
% load the model data before running
close all
clear all 

ref_file = 2;

% C_nuc_rate
data_files = [27:31];

% % C_qdot_lw
% data_files = [18:21];
% 
% % C_coalescence (collision eff. only)
% % add in 9 to give one where C = 0
% data_files = [9,22:25];




cd('/Users/jez/Desktop/stuff from nitro black/model_results');

data = load(['bubble_sim_data' num2str(ref_file)]);

P_ref = data.P;
t_ref = data.t;
for i = 1:length(data_files)
    clear data
    filename = ['bubble_sim_data' num2str(data_files(i))];
    data = load(filename);
    
        P_compare = interp1(data.t, data.P, t_ref);

    
    figure(1)
    hold on
    plot(data.t, data.P/1e6, 'k')
    
        figure(2)
    hold on
    plot(t_ref, (P_compare - P_ref), 'k')

    
end

figure(1)
xlabel('Time [s]')
ylabel('Pressure [MPa]')
make_text_big('shift_start')
make_text_big

figure(2)
xlabel('Time [s]')
ylabel('Pressure Difference [Pa]')
make_text_big('shift_start')
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