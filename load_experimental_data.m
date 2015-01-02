function P_out = load_experimental_data(t_sim)

load data_for_finding_constants

P_out = interp1(t/t(end), P, t_sim/t_sim(end));