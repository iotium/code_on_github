function node_level = get_node_levels(V_l, V_bubi, V_node, node_level_guess)
IC = sum(node_level_guess);

x = fzero(@(x) error_fn(x), IC);

node_level = convert_x_to_NL(x);


    function E = error_fn(x)
        
        NL = convert_x_to_NL(x);
        
        if sum(size(NL) == size(V_bubi)) ~= 2
            disp('uh oh, tank is overfilled in node level calculation')
        end
        
        E = (V_l/V_node - sum(NL) + sum(NL.*V_bubi));
        
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