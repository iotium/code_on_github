function varargout = unpack_y(varargin)

y = varargin{1};
constants = varargin{2};

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

if nargin == 2
    
    variables.m_tg = m_tg;
    variables.U_tg = U_tg;
    variables.T_gw = T_gw;
    variables.m_l = m_l;
    variables.T_l = T_l;
    variables.T_lw = T_lw;
    
    variables.g_q = g_q;
    variables.w_q = w_q;
    
    varargout{1} = variables;
    
else
    
    ind_max_rel_err = varargin{3};
   rel_err = varargin{4};
    
    fprintf('max error is at ')
    if ind_max_rel_err < 4+2*N_rw
        % it's a wall temp or a liquid/ullage property
        if ind_max_rel_err == 1
            % it's ullage mass
            fprintf('ullage mass')
        elseif ind_max_rel_err == 2
            % it's ullage internal energy
            fprintf('ullage U')
        elseif ind_max_rel_err >=3 && ind_max_rel_err <= 2+N_rw
            % it's an ullage wall temp
            fprintf('ullage wall temp')
        elseif ind_max_rel_err == 3+N_rw
            % it's liquid mass
            fprintf('liquid mass')
        elseif ind_max_rel_err == 4+N_rw
            % it's liquid temp
            fprintf('liquid temp')
        else
            % it's liquid wall temp
            fprintf('liquid wall temp')
        end
    else
        % it's a weight or abscissa
        
        if sum(ind_max_rel_err == i_w) == 1
            % it's a weight
            
            % the index of the weight
            ind_w_err = ind_max_rel_err - (4+2*N_rw);
            
            [ind_node, ind_weight] = ind2sub([ N_nodes, N_ab], ind_w_err);
            
            fprintf('node %0.d, weight no. %0.d', ind_node, ind_weight)
        else
            % it's an abscissa
            
            % the index of the abscissa
            ind_g_err = ind_max_rel_err - i_w(end);
            
            [ind_node, ind_abs] = ind2sub([ N_nodes, N_ab], ind_g_err);
            
            fprintf('node %0.d, weighted abscissa no. %0.d', ind_node, ind_abs)
        end
    end
    fprintf(' rel error is %0.3g, abs(y) = %0.3g\n', rel_err, abs(y(ind_max_rel_err)) )
end


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