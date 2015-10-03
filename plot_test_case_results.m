% plotting data for thesis
% load the model data before running
close all

[P_exp, T_lw_out_exp, LL_exp] = load_experimental_data(t, specified_case);

m_out = cumtrapz(t, mdot_out_liq + mdot_out_vap);
mh_out = cumtrapz(t, (mdot_out_liq.*h_l + mdot_out_vap.*h_tg_sat));

T_l = y(4+N_rw,:);
m_l = y(3+N_rw,:);

m_tg = y(1,:);

T_lw_in = y(5+N_rw,:);
T_gw_in = y(3,:);
T_lw_out = y(4+2*N_rw,:);

% P vs t
figure(1)
plot(t,P/1e6,'k', t, P_exp/1e6,'b')
legend('Model','Experiment')
ylabel('Pressure [MPa]')
xlabel('Time [s]')
make_text_big('shift_start')


% fill level vs t
figure(2)
plot(t, fill_level*100, 'k', t, LL_exp*100, 'b')
xlabel('Time [s]')
ylabel('Fill Level [%]')
legend('Model','Experiment')
make_text_big

% % find P at 25% fill level
% t_LL_25 = interp1(fill_level, t, 0.25);
% P_LL_25 = interp1(fill_level, P, 0.25);

figure(3)
plot(t,T_l-273.15,'k',t,T_s-273.15,'k--')
legend('Liquid','T_{sat}(P) = T_{tg}')
ylabel('Temperature [C]')
xlabel('Time [s]')
make_text_big

figure(4)
plot(t, mdot_out_liq + mdot_out_vap, 'k')
hold on
plot(t, gas_holdup_injector, 'b')
xlabel('Time [s]')
ylabel('Mass Flow Rate [kg/s]')
legend('Mass Flow Rate','Vapor Mass Fraction')
make_text_big

figure(5)
plot(t, x_tg,'k')
xlabel('Time [s]')
ylabel('Ullage Vapor Mass Fraction')
make_text_big('shift_start')


t_plot = linspace(0,1,10)*t(end);


% find sauter mean diameter, void fraction, and number density
% versus space, at time indices of t_plot
SMD = NaN*ones(N_nodes,length(t_plot));
alpha = SMD;
ND = SMD;
% dx_nodes = 1/(N_nodes-1);
x_nodes = linspace(0,1,N_nodes)';
r_nodes = linspace(0,1,N_rw);


for i = 1:length(t_plot)
    [~,i_plot] = min(abs(t - t_plot(i)));
    ind_full = 1:(N_full(i_plot));
    SMD(ind_full,i) = mom(ind_full,4,i_plot)./mom(ind_full,3,i_plot);
    alpha(ind_full,i) = V_bubi(ind_full,i_plot);
    ND(ind_full,i) = mom(ind_full,1,i_plot);
    
    figure(6)
    hold on
    plot(x_nodes, SMD(:,i), 'k')

    figure(7)
    hold on
    plot(x_nodes, alpha(:,i), 'k')

    figure(8)
    hold on
    plot(x_nodes, ND(:,i), 'k')
    
    figure(9)
    hold on
    plot(r_nodes, T_lw(:,i_plot)-273.15, 'k')
    
    figure(10)
    hold on
    plot(r_nodes, T_gw(:,i_plot)-273.15, 'k')

end

figure(6)
xlabel('Normalized Height')
ylabel('Sauter Mean Diameter [m]')
make_text_big

figure(7)
xlabel('Normalized Height')
ylabel('Void Fraction')
make_text_big

figure(8)
xlabel('Normalized Height')
ylabel('Bubble Number Density [#/m^3]')
make_text_big

figure(9)
xlabel('Normalized Radius')
ylabel('Temperature [C]')
make_text_big

figure(10)
xlabel('Normalized Radius')
ylabel('Temperature [C]')
make_text_big

figure(11)
plot(t, gas_holdup_injector, 'k')
xlabel('Time [s]')
ylabel('Outlet Vapor Mass Fraction')






% 
% figure(6)
% plot(t, gas_holdup_injector, 'k')
% xlabel('Time [s]')
% ylabel('Outlet Vapor Mass Fraction')
% make_text_big


% figure(4)
% plot(t,T_lw_in,'k-',t,T_gw_in,'b--',t,T_lw_out,'r-.',t,T_lw_out_exp+273.15,'k:')
% title('wall temp')
% xlabel('Time [s]')
% legend('liquid','vapor','outside liquid','experimental outside liquid')
% title('wall temp')