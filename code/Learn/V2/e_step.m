function q=e_step(A,hidden_indexes,hidden_vals,theta_c,training_data,p_x_observed)
data_size=size(training_data,1);
hidden_size=size(hidden_vals,1);
% q hols the p(x"miss" | x"observed") for all training example i and every
% possible value of x"miss"
q=zeros(data_size,hidden_size);
for i=1:data_size
    for j=1:hidden_size
        training_data(i,hidden_indexes)=hidden_vals(j,:);
        q(i,j)=calc_p_all_x(A,training_data(i,:),theta_c)/p_x_observed(i,1);
    end
end
end