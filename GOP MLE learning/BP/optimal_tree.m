function [ll, time] =  optimal_tree(adj_mat, nodes_values, training_data, opt)
% opt: is structure holding all the options 
% opt.tol: is the tolerance
% opt.max_ite: is the maximum no. of iterations
% opt.LR: is the learning rate
% nodes_values: is a vector where nodes_values(i) = the max value node i
% can get

tic

iterations = 1;
ll_new = -inf;
ll_old = 0;
nodes_count = size(adj_mat,1);
training_size = size(training_data, 1);
% intialize theta and build the tree
[ tree, theta, lambda] = create_tree(adj_mat, nodes_values);
% get the count of hidden and observed valus
observed_variables = size(training_data, 2);
hidden_variables = nodes_count - observed_variables;

while abs(ll_new - ll_old) < opt.tol || iterations < opt.max_ite
    ll_old = ll_new;
    ll_new = 0;
    iterations = iterations + 1;
    % given current theta calculate z and beliefs when all nodes are not
    % clamped
    [z_beliefs, z] = inference(tree,theta,ones(1, nodes_count) * -1, 0, theta);
    beliefs = z_beliefs;
    for i = 1:training_size
        current_ll = 0;
        fixed_nodes = ones(1, nodes_count) * -1;
        fixed_nodes(1, 1:observed_variables) = training_data(i,:); 
        if i == 1
            [beliefs, current_ll] = inference(tree,theta,fixed_nodes, 0, beliefs);
        else
            [beliefs, current_ll] = inference(tree,theta,fixed_nodes, 1, beliefs);
        end
        ll_new = ll_new + log(current_ll/z);
    end
    ll_new = ll_new / training_size;
    theta = update_theta(theta, z_beliefs, beliefs, lambda, opt.LR, training_size);
    fprintf('ll = %f   , iterations = %i \n', ll_new, iterations);
end
ll = ll_new;
time = toc;
end