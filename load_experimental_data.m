function [P_out, T_surface_out] = load_experimental_data(t_sim)

load their_data

% P = N2O_11.p;
% t = N2O_11.t;
% 
% P = N2O_2.p;
% t = N2O_2.t;

% P = quartz24.p;
% t = quartz24.t;

% P = quartz48.p;
% t = quartz48.t;

% P = glass_gauge236.p;
% t = glass_gauge236.t;

% P = glass_gauge291.p;
% t = glass_gauge291.t;
% T_surface = glass_gauge291.T_surface;

P = glass_gauge289.p;
t = glass_gauge289.t;
T_surface = glass_gauge289.T_surface;

% P_out = interp1(t/t(end), P, t_sim/t_sim(end));
P_out = interp1(t, P, t_sim);
T_surface_out = interp1(t, T_surface, t_sim);