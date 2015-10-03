function out = sum_over_nodes(prop_in, node_level, V_node)
% calculates the net amount of prop_in in all the nodes.
% prop_in is something/(unit volume).
% node_level is the fill level of each node (0-1)
% adaptive mesh refinement version

if length(prop_in) ~= length(node_level)
    error('input and node levels are different lengths')
end

out = sum( node_level .* V_node .* prop_in);


%{
% V_node is the nominal node volume (top and bottom nodes V = 1/2 V_node)

if length(prop_in) ~= length(node_level)
    error('input and node levels are different lengths')
end

out = V_node/2*prop_in(1)*node_level(1) ...
    + sum(node_level(2:end-1).*prop_in(2:end-1)*V_node) ...
    + V_node/2*prop_in(end)*node_level(end);
%}