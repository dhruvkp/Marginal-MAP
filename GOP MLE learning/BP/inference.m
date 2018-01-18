function [beliefs, result]= inference(input_tree,theta,fixed_nodes, opt, beliefs)
% opt: is option takes values 0 or 1, where zeros means calc only current
% beliefs, and 1 means add to the previous ones to
% fixed node is the set of predefined node for example to calculate 
% Z all entries will be set to -1

% calc bottom up msgs
input_tree = calc_up_msgs(input_tree,theta,fixed_nodes);
% calc up down msgs
input_tree = calc_down_msgs(input_tree,theta,fixed_nodes);
% calc the product of all incomming messages to all nodes
products = calc_msg_prod(input_tree);
result=0;
% now all the messages have been delivered to the root
% summing them up will give us the required probability
for j=1:input_tree{1,3}
    multi=1;
    for i=1:size(input_tree{1,2},2)
        multi=multi*input_tree{1,2}(j+1,i);
    end
    result=result+multi;
end

% find the edges in the tree
[i, j] = find(~cellfun(@isempty,beliefs));
indexes = [i, j];

% calculate beliefs on every edge  
for k = 1:size(indexes, 1)
    % calc the product of all incomming messages to the two nodes exept the
    % one comming from the other
    current_parent = indexes(k,1);
    current_child = indexes(k,2);
    parent_prod = products{current_parent,1} ./ ...
        input_tree{current_parent,2}(2:size(input_tree{current_parent,2},1), find(input_tree{current_parent,2}(1,:)==current_child));
    child_prod = products{current_child,1} ./ ...
        input_tree{current_child,1}(2:size(input_tree{current_child,1},1), find(input_tree{current_child,1}(1,:)==current_parent));
    parent_prod(isnan(parent_prod)) = 0;
    child_prod(isnan(child_prod)) = 0;
    if opt == 0
        beliefs{current_parent, current_child} = parent_prod * child_prod';
    else
        beliefs{current_parent, current_child} = beliefs{current_parent, current_child} + parent_prod' * child_prod;
    end
    
    beliefs{current_parent, current_child} = beliefs{current_parent, current_child} ...
                                            .* theta{current_parent, current_child};
end
end