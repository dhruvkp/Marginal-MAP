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
count=0;
while ~isequal(theta_c,temp_theta)
    count=count+1;
    temp_theta=theta_c;
    q=e_step(A,hidden_indexes,hidden_vals,theta_c,training_data,p_x_observed);
    q=q ./ calc_normal(A,binary_combinations(number_variables),theta_c);
    theta_c=m_step(A,hidden_indexes,hidden_vals,theta_shape,training_data,q);
end
end