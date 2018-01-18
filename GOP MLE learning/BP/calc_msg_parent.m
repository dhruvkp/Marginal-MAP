function [node_list input_tree]=calc_msg_parent(input_tree, node_no, theta, fixed_nodes)
% get parents of the node
parent_list=input_tree{node_no,1}(1,:);
node_list=[];
% for all the parents do calculate the message
for j=1:size(parent_list,2)
    current_parent=parent_list(j);
    for xj=1:input_tree{current_parent,3}
        % if the value of the parent is predefined set all the messages
        % about other values to zero
        if fixed_nodes(1,current_parent)~= -1 
            if fixed_nodes(1,current_parent)~= xj
                input_tree{current_parent,2}(xj+1,find(input_tree{current_parent,2}(1,:)==node_no))=0;
               continue;
            end
        end
        sum=0;
        % for all possible values of the child calculate the part of the
        % message regard this value
        for xi=1:input_tree{node_no,3}
            message_mult=1;
            % if the child value is predefined the message regard all other
            % values will be 1
            if fixed_nodes(1,node_no)~= -1 
                if fixed_nodes(1,node_no)~= xi
                   continue;
                end
            end
            % message regard any value is the multible of all the incoming
            % messages about this value and the beliefs of the parent,
            % child values
            for child=1:size(input_tree{node_no,2},2)
                message_mult=message_mult* input_tree{node_no,2}(xi+1,child);
            end
            % finally the message to the parent for a specific value is the over all 
            %the incomming messages regard this value
            sum=sum+(exp(theta{current_parent,node_no}(xj,xi))*message_mult);
        end
        input_tree{current_parent,2}(xj+1,find(input_tree{current_parent,2}(1,:)==node_no))=sum;
      
    end
    node_list=[node_list current_parent];
end
end