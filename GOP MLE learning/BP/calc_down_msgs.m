function  input_tree = calc_down_msgs(input_tree,theta,fixed_nodes)
% fixed node is the set of predefined node for example to calculate 
% Z all entries will be set to -1

% starting from the leafs up
leaf_list=get_leafs(input_tree);
node_list=[1];
while size(node_list,2)>=1
    current_node=node_list(1);
    % get the first element in the list
    node_list=node_list(:,2:size(node_list,2));
    if ismember(current_node, leaf_list)
        continue;
    end
    % calculate the node message to it's parents
    [result input_tree]=calc_msg_child1(input_tree, current_node, theta,fixed_nodes);
    % add the parents to our list
    node_list=[node_list result];
    % delete redundant
    node_list=unique(node_list);
    %prevent stuking in the root
    node_list(ismember(node_list,1))=[];
end

end