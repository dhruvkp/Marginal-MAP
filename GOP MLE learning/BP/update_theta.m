function theta = update_theta(theta, z_beliefs, beliefs, lambda, LR, training_size)

% find the edges in the tree
[i, j] = find(~cellfun(@isempty,beliefs));
indexes = [i, j];

% update theta on every edge  
for k = 1:size(indexes, 1)
    % calc the product of all incomming messages to the two nodes exept the
    % one comming from the other
    current_parent = indexes(k,1);
    current_child = indexes(k,2);
    theta{current_parent, current_child} = theta{current_parent, current_child} + ...
                                           LR * (beliefs{current_parent, current_child} - ...
                                           (training_size * z_beliefs{current_parent, current_child}) - ...
                                           (lambda * theta{current_parent, current_child}));
end
end