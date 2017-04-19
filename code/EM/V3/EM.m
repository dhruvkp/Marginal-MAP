function [theta_c] = EM( A,input_indexes,output_indexes,hidden_indexes,x_domain,training_data)
number_variables=size(A,1);
domain_size=size(x_domain,2);
theta_c=rand(number_variables,number_variables,domain_size,domain_size);
for i=1:number_variables
    for j=1:number_variables
        if A(i,j)>0
            for x_i=x_domain
                for x_j=x_domain
                    theta_c(i,j,x_i,x_j)=(x_i+x_j)/12;
                end
            end
        end
    end
end
% for i=1:domain_size
%     for j=1:domain_size
%         theta_c(:,:,i,j)=theta_c(:,:,i,j).* A;
%     end
% end

theta_shape=size(theta_c);
hidden_vals=binary_combinations(size(hidden_indexes,1));
p_x_observed=calc_p_x_observed(A,hidden_indexes,hidden_vals,theta_c,training_data);
temp_theta=zeros(theta_shape);
% all posible values for the hidden variables compined
h_y_combination=binary_combinations(size(hidden_indexes,1)+size(output_indexes,1));
h_y_vals=zeros(size(h_y_combination,1),size(A,1));
h_y_vals(:,hidden_indexes)=h_y_combination(:,1:size(hidden_indexes,1));
h_y_vals(:,output_indexes)=h_y_combination(:,1:size(output_indexes,1));
[training_unique,~,u_id] = unique(training_data(:,input_indexes), 'rows'); % // find unique rows and their unique id
occurrences = histc(u_id, unique(u_id)); % // count occurrences of unique ids
count=0;
t=0;
while t<100 %~isequal(theta_c,temp_theta)
    count=count+1;
    temp_theta=theta_c;
    q=e_step(A,hidden_indexes,hidden_vals,theta_c,training_data,p_x_observed);
    q=q ./ calc_normal(A,binary_combinations(number_variables),theta_c);
    theta=m_step(A,hidden_indexes,h_y_vals,training_unique,occurrences,input_indexes,hidden_vals,theta_shape,training_data,q);
    y1 = isnan(theta);
    if size(theta,1)==0 || any(isnan(theta(:))) 
        break;
    else
    theta_c=theta;
    end
    tt=theta_c-temp_theta;
    t=t+1;
end
end