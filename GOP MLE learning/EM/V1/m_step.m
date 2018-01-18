function theta=m_step(A,hidden_indexes,hidden_vals,theta_shape,training_data,q)
theta=zeros(theta_shape);
number_variables=size(A,1);
data_size=size(training_data,1);
hidden_size=size(hidden_vals,1);
for i=1:data_size
    for j=1:hidden_size
        training_data(i,hidden_indexes)=hidden_vals(j,:);
        for k=1:number_variables
            for l=1:number_variables
                if A(k,l)>0 
                    theta(k,l,training_data(i,k),training_data(i,l))= ...
                        theta(k,l,training_data(i,k),training_data(i,l))+q(i,j);
                        
               end
            end
        end
    end
end
end