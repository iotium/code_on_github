function variables = unpack_y(y, constants)

N_rw = constants.N_rw;

N_ab = constants.N_ab;

N_mom = N_ab*2;

N_nodes = constants.N_nodes;

% indices of the weights
i_w = (4 + 2*N_rw) + [1:(N_nodes*N_ab)];

w_q = y(i_w);

% reshape so each ROW is one node
w_q = reshape(w_q, N_ab, N_nodes)';

% indices of the abscissas
i_g = i_w(end) + [1:(N_nodes*N_ab)];

g_q = y(i_g);

% reshape so each ROW is one node
g_q = reshape(g_q, N_ab, N_nodes)';

m_tg = y(1);
U_tg = y(2);
T_gw = y(3 : 2+N_rw);
m_l = y(3+N_rw);
T_l = y(4+N_rw);
T_lw = y(5+N_rw : 4+2*N_rw);

variables.m_tg = m_tg;
variables.U_tg = U_tg;
variables.T_gw = T_gw;
variables.m_l = m_l;
variables.T_l = T_l;
variables.T_lw = T_lw;

variables.g_q = g_q;
variables.w_q = w_q;



% # of variables per node: N_mom (weights + weighted abscissas)
% total number of variables: 6 + N_mom*N_nodes

% 1st node is at bottem of tank
% weights and abscissas are sorted by #, not by node
% so first all the 1st abscissas are given for each node, then the 2nd...
% IE like this:

% 1st ab for 1st node
% 1st ab for 2nd node
% 1st ab for 3rd node
% 2nd ab for 1st node
% 2nd ab for 2nd node
% ...

% scratch that I'll do the opposite
% 1st weight for 1st node
% 2nd weight for 1st node
% 1st weight for 2nd node
% 2nd weight for 2nd node
% ...
% 1st ab. for 1st node

% N_nodes = (length(y) - (4 + 2*N_rw)) / N_mom;
% 
% % indices of the weights
% i_w = (4 + 2*N_rw) + [1:(N_nodes*N_ab)];
% 
% w_q = y(i_w);
% 
% % reshape so each ROW is one node
% w_q = reshape(w_q, N_ab, N_nodes)';
% 
% % indices of the abscissas
% i_g = i_w(end) + [1:(N_nodes*N_ab)];
% 
% g_q = y(i_g);
% 
% % reshape so each ROW is one node
% g_q = reshape(g_q, N_ab, N_nodes)';
% 
% m_tg = y(1);
% U_tg = y(2);
% T_gw = y(3 : 2+N_rw);
% m_l = y(3+N_rw);
% T_l = y(4+N_rw);
% T_lw = y(5+N_rw : 4+2*N_rw);