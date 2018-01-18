function [node_list input_tree]=calc_msg_child(input_tree, node_no, theta, fixed_nodes)
% get parents of the node
childs_list=input_tree{node_no,2}(1,:);
node_list=[];
% for all the childs do calculate the message
for j=1:size(childs_list,2)
    current_child=childs_list(j);
    for xj=1:input_tree{current_child,3}
        % if the value of the child is predefined set all the messages
        % about other values to zero
        if fixed_nodes(1,current_child)~= -1 
            if fixed_nodes(1,current_child)~= xj
                input_tree{current_child,1}(xj+1,find(input_tree{current_child,1}(1,:)==node_no))=0;
               continue;
            end
        end
        sum=0;
        % for all possible values of the child calculate the part of the
        % message regard this value
        for xi=1:input_tree{node_no,3}
            message_mult=1;
            if node_no == 1
               message_mult = prod(input_tree{1,2}(xi+1, :)) / ...
                              input_tree{1,2}(xi+1, find(input_tree{1,2}(1,:)==current_child));
            end
            % if the child value is predefined the message regard all other
            % values will be 1
            if fixed_nodes(1,node_no)~= -1 
                if fixed_nodes(1,node_no)~= xi
                   continue;
                end
            end
            % message regard any value is the multible of all the incoming
            % messages about this value and the beliefs of the child,
            % child values
            for child=1:size(input_tree{node_no,1},2)
                message_mult=message_mult* input_tree{node_no,2}(xi+1,child);
            end
            % finally the message to the child for a specific value is the over all 
            %the incomming messages regard this value
            sum=sum+(theta{node_no, current_child}(xi,xj)*message_mult);
        end
        input_tree{current_child,1}(xj+1,find(input_tree{current_child,1}(1,:)==node_no))=sum;
    end
    node_list=[node_list current_child];
end
end