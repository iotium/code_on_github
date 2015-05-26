function node_level = get_node_levels(V_l, V_bubi, V_node, node_level_guess)
% gets the node levels based on the current liquid and bubble volumes
% this is the version for points on the boundaries (ie FD scheme, not FV)

% IC = sum(node_level_guess);
% 
% x = fzero(@(x) error_fn(x), IC);
% % x = fsolve(@(x) error_fn(x), IC, optimset('Jacobian','on','display','off'));
% 
% x_NI = non_iterative_version;
% 
% if abs(x_NI - x) > 1e-6
%     disp('non-iterative version didn''t work')
% end

x = non_iterative_version;

node_level = convert_x_to_NL(x);

    function x = non_iterative_version
        % for a description, see page 101 of book #6
        
        N = length(V_bubi);
        
        bub_term = cumsum(V_bubi(2:N-1));
        bub_term = [0; bub_term(:)];
        
        N_full_term = [1:N-2]';
        
        % for N_full = 0
        x_top(1) = V_l*2/(V_node * (1 - V_bubi(1)) );
        
        % for N_full > 0 and < N-1 
        x_top(2:N-1) = (V_l/V_node + V_bubi(1)/2 + bub_term(1:N-2) - N_full_term + 0.5)./(1 - V_bubi(2:N-1));
        
        % for N_full = N-1
        x_top(N) =   2*(V_l/V_node + V_bubi(1)/2 + bub_term(N-1) - (N-1) + 0.5)./(1 - V_bubi(N));
        
%         x_top = (V_l/V_node - N_full_term + bub_term)./(1 - V_bubi);

        ind_correct = find( x_top > 0 & x_top < 1);
        N_full_term = [0:N-1];
        N_full = N_full_term(ind_correct);
        
        if length(N_full) > 1
            error('node level problem')
        end
        
        x_top = x_top(ind_correct);
        
        x = x_top + N_full;
    end


    function [E] = error_fn(x)
        % used to get the node levels in an iterative way
        % the "non iterative version" is much faster
        
        
        % x is the actual number (eg 7.65)
        
        NL = convert_x_to_NL(x);
        
        if sum(size(NL) == size(V_bubi)) ~= 2
            disp('uh oh, tank is overfilled in node level calculation')
        end
        
        E = (V_l/V_node - x + sum(NL.*V_bubi));
        
        % jacobian
        J = -1 + V_bubi(floor(x) + 1);
        
    end

    function NL = convert_x_to_NL(x)
        % convert from a single number (eg 7.65) to a value for each node:
        % (eg [1 1 1 1 1 1 1 0.65 0 0 0 ...])
        
        N_nodes = length(V_bubi);
        
        NL = zeros( N_nodes, 1);
        
        if x > 1
            
            NL(1:floor(x)) = 1;
            
        end
        
        if ceil(x) <= 0
            disp('problem in getting node levels - it went negative')
        end
        
        if x > length(V_bubi)
           
            disp('problem in getting node levels - fill level got larger than tank volume')
            NL = ones(size(V_bubi));
        else
            
            NL(ceil(x)) = rem(x,1);
        end
        
    end

end