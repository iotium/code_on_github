function [P_out, T_surface_out, LL_out] = load_experimental_data(t_sim, specified_case)

load their_data

switch specified_case
    case 5
        P = N2O_2.p;
        t = N2O_2.t;
    case 6
        P = N2O_11.p;
        t = N2O_11.t;
        
    case 7
        
        P = quartz24.p;
        t = quartz24.t;
    case 8
        P = quartz48.p;
        t = quartz48.t;
        t_LL =quartz48.t_LL;
        LL = quartz48.LL;
    case 9
        P = glass_gauge236.p;
        t = glass_gauge236.t;
        t_LL = glass_gauge236.t_LL;
        LL = glass_gauge236.LL;
    case 10
        P = glass_gauge291.p;
        t = glass_gauge291.t;
        T_surface = glass_gauge291.T_surface;
        t_LL = glass_gauge291.t_LL;
        LL = glass_gauge291.LL;
        
    case 11
        P = glass_gauge289.p;
        t = glass_gauge289.t;
        T_surface = glass_gauge289.T_surface;
        t_LL = glass_gauge289.t_LL;
        LL = glass_gauge289.LL;
    case 12
        P = glass_gauge257.P;
        t = glass_gauge257.t;
        t_LL = glass_gauge257.t_LL;
        LL = glass_gauge257.LL;
    otherwise
                P = glass_gauge289.p;
        t = glass_gauge289.t;
        T_surface = glass_gauge289.T_surface;
        t_LL = glass_gauge289.t_LL;
        LL = glass_gauge289.LL;
end

% P_out = interp1(t/t(end), P, t_sim/t_sim(end));
P_out = interp1(t, P, t_sim);

if exist('T_surface','var')
    T_surface_out = interp1(t, T_surface, t_sim);
else
    T_surface_out = nan*ones(size(t_sim));
end

if exist('LL','var')
    LL_out = interp1(t_LL, LL, t_sim);
else
    LL_out = nan*ones(size(t_sim));
end