function [node_list input_tree]=calc_msg_parent1(input_tree, node_no, theta, fixed_nodes)
% get parents of the node
parent_list=input_tree{node_no,1}(1,:);
node_list=[];
% for all the parents do calculate the message
child_values = zeros(1, input_tree{node_no,3})';
if fixed_nodes(1,node_no)~= -1
    child_values( fixed_nodes(1,node_no), 1) = 1;
else
    child_values = ones(1, input_tree{node_no,3})';
end
for j=1:size(parent_list,2)
    current_parent=parent_list(j);
    parent_values = zeros(1, input_tree{current_parent,3})';
    if fixed_nodes(1,current_parent)~= -1
        parent_values( fixed_nodes(1,current_parent),1) = 1;
    else
        parent_values = ones(1, input_tree{current_parent,3})';
    end
    if isempty(input_tree{node_no,2}(2:size(input_tree{node_no,2},1),:))
        messages = (exp(theta{current_parent,node_no} )* child_values) .* parent_values;
    else
        messages = (exp(theta{current_parent,node_no} )* (child_values .* prod(input_tree{node_no,2}(2:size(input_tree{node_no,2},1),:), 2))) .* parent_values;
    end
    input_tree{current_parent,2}(2:size(input_tree{current_parent,2},1),find(input_tree{current_parent,2}(1,:)==node_no))=messages;
      
    
    node_list=[node_list current_parent];
end
end