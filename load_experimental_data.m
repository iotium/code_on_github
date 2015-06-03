function P_out = load_experimental_data(t_sim)

load their_data

% P = N2O_11.p;
% t = N2O_11.t;

P = N2O_2.p;
t = N2O_2.t;

% P_out = interp1(t/t(end), P, t_sim/t_sim(end));
P_out = interp1(t, P, t_sim);