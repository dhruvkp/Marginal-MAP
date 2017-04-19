function theta_c=initialize_theta(A,number_variables,x_domain)
domain_size=size(x_domain,1);
theta_c=zeros(number_variables,number_variables,domain_size,domain_size);
for i=1:number_variables
    for j=1:number_variables
        if A(i,j)>0
            for x_i=x_domain
                for x_j=x_domain
                    if j<i
                        theta_c(i,j,x_i,x_j)=theta_c(j,i,x_i,x_j);
                    else
                        theta_c(i,j,x_i,x_j)=rand();
                    end
                end
            end
        end
    end
end
end