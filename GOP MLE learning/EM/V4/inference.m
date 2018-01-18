function result=inference(input_tree,theta,fixed_nodes)
% fixed node is the set of predefined node for example to calculate 
% Z all entries will be set to -1

% starting from the leafs up
leaf_list=get_leafs(input_tree);
node_list=leaf_list;
while size(node_list,2)>=1
    current_node=node_list(1);
    % get the first element in the list
    node_list=node_list(:,2:size(node_list,2));
    % calculate the node message to it's parents
    [result input_tree]=calc_msg_parent(input_tree, current_node, theta,fixed_nodes);
    % add the parents to our list
    node_list=[node_list result];
    % delete redundant
    node_list=unique(node_list);
    %prevent stuking in the root
    node_list(ismember(node_list,1))=[];
end
result=0;
% now all the messages have been delivered to the root
% summing them up will give us the required probability
for j=1:2
    multi=1;
    for i=1:size(input_tree{1,2},2)
        multi=multi*input_tree{1,2}(j+1,i);
    end
    result=result+multi;
end
end