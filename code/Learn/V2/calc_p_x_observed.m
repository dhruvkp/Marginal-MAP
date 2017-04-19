% calculate the marginalized probability of the observed x of every training example 
function p_list=calc_p_x_observed(A,hidden_indexes,hidden_vals,theta_c,training_data)
% input: adjacency matrix, clique parameters, burn-in iterations, number of
% samples required
% output: samples
data_size=size(training_data,1);
p_list=zeros(data_size,1);
for i=1:data_size
    for j=1:size(hidden_vals,1)
        training_data(i,hidden_indexes)=hidden_vals(j,:);
        p_list(i,1)=p_list(i,1)+calc_p_all_x(A,training_data(i,:),theta_c);
    end
end
end