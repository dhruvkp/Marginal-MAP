function [node_list input_tree]=calc_msg_child1(input_tree, node_no, theta, fixed_nodes)
% get parents of the node
childs_list=input_tree{node_no,2}(1,:);
node_list=[];
% for all the parents do calculate the message
for j=1:size(childs_list,2)
    current_child=childs_list(j);
    child_values = zeros(1, input_tree{current_child,3})';
    if fixed_nodes(1,current_child)~= -1
        child_values(fixed_nodes(1,current_child), 1) = 1;
    else
        child_values = ones(1, input_tree{current_child,3})';
    end
    product = 1;
    if node_no == 1
        product = prod(input_tree{node_no,2}(2:size(input_tree{node_no,2},1),:),2)./ ...
                          input_tree{1,2}(2:size(input_tree{1,2},1), find(input_tree{1,2}(1,:)==current_child));
    else
        product = prod(input_tree{node_no,2}(2:size(input_tree{node_no,2},1),:),2);
    end
    product(isnan(product)) = 0;
    messages = (theta{node_no, current_child}' *  product) .* child_values;
    input_tree{current_child,1}(2:size(input_tree{current_child,1},1),find(input_tree{current_child,1}(1,:)==node_no))=messages;
      
    
    node_list=[node_list current_child];
end
end