function result_theta=m_step(A,theta,h_y_vals,training_unique,occurrences,input_indexes)
result_theta=zeros(size(theta));
number_variables=size(A,1);
data=h_y_vals;
data_size=size(data,1);
training_size=size(training_unique,1);
% calculating theta= argmax(theta) sum_over_data"i"(sum_over_all_possible_hidden"j"(q(i,j)* log(p(x"observed", j) )))
temp_theta=ones(size(theta));
for f=1:training_size
    z=0;temp_theta=ones(size(theta));
    for i=1:data_size
        data(i,input_indexes)=training_unique(f,:);
        current=calc_p_all_x(A,data(i,:),theta);
        if current == Inf
            display(i);
        end
        z=z+current;
            for k=1:number_variables
                for l=k:number_variables

                    if A(k,l)>0 
                        % close form to update theta
                        result_theta(k,l,data(i,k),data(i,l))= ...
                            result_theta(k,l,data(i,k),data(i,l))+(current*occurrences(f,1));
                        temp_theta(k,l,data(i,k),data(i,l))=0;
                        temp_theta(l,k,data(i,l),data(i,k))=0;
                        result_theta(l,k,data(i,l),data(i,k))= ...
                            result_theta(k,l,data(i,k),data(i,l));
%                         result_theta(k,l,data(i,l),data(i,k)) = ...
%                             result_theta(l,k,data(i,k),data(i,l));
                   end
                end

             end
    end
    temp_theta(find(temp_theta==0))=z;
    result_theta=result_theta./temp_theta;
end
result_theta=result_theta./sum(occurrences);
%result_theta=log(result_theta);
end