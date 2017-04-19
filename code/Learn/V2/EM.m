function [theta_c] = EM( A,input_indexes,output_indexes,hidden_indexes,x_domain,training_data)
number_variables=size(A,1);
% initialize theta
theta_c=initialize_theta(A,number_variables,x_domain);

theta_shape=size(theta_c);
% all posible values for the hidden variables compined
h_y_combination=binary_combinations(size(hidden_indexes,1)+size(output_indexes,1));
h_y_vals=zeros(size(h_y_combination,1),size(A,1));
h_y_vals(:,hidden_indexes)=h_y_combination(:,1:size(hidden_indexes,1));
h_y_vals(:,output_indexes)=h_y_combination(:,1:size(output_indexes,1));
[training_unique,~,u_id] = unique(training_data(:,input_indexes), 'rows'); % // find unique rows and their unique id
occurrences = histc(u_id, unique(u_id)); % // count occurrences of unique ids
% for each training case calculate the probability of it's observed part
% % to be used in calculating q in every iteration
% p_x_observed=calc_p_x_observed(A,hidden_indexes,hidden_vals,theta_c,training_data);

% temp theta to cheack convergence
temp_theta=zeros(theta_shape); 
tt=theta_c-temp_theta;
% while not converge
while sum(sum(sum(sum(abs(tt)))))>1*10^-10
    temp_theta=theta_c;
    % E step
%     q=e_step(A,hidden_indexes,hidden_vals,theta_c,training_data,p_x_observed);
%     % renormalize q
%     q=q ./ calc_normal(A,binary_combinations(number_variables),theta_c);
    %M step
   theta_c=m_step(A,theta_c,h_y_vals,training_unique,occurrences,input_indexes);
   tt=theta_c-temp_theta;
   %tt
end
end