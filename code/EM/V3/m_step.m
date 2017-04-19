function theta=m_step(A,hidden_indexes,h_y_vals,training_unique,occurrences,input_indexes,hidden_vals,theta_shape,training_data,q)
f=zeros(theta_shape);
theta=zeros(theta_shape);
theta1=zeros(theta_shape);
number_variables=size(A,1);
data_size=size(training_data,1);
hidden_size=size(hidden_vals,1);
for i=1:data_size
    for j=1:hidden_size
        training_data(i,hidden_indexes)=hidden_vals(j,:);
        for k=1:number_variables
            for l=1:number_variables
                if A(k,l)>0 
                    f(k,l,training_data(i,k),training_data(i,l))= ...
                        f(k,l,training_data(i,k),training_data(i,l))+q(i,j);
                        
               end
            end
        end
    end
end
data=h_y_vals;
data_size=size(data,1);
training_size=size(training_unique,1);
for n=1:training_size
    for j=1:hidden_size
        for i=1:data_size
            data(i,input_indexes)=training_unique(n,:);
            current=calc_p_all_x(A,data(i,:),theta);
            if current == Inf
                display(i);
            end
                for k=1:number_variables
                    for l=k:number_variables

                        if A(k,l)>0 
                            % close form to update theta
                            f(k,l,data(i,k),data(i,l))= ...
                                f(k,l,data(i,k),data(i,l))+(current*occurrences(n,1)*q(n,j));
                            f(l,k,data(i,l),data(i,k))= ...
                                f(k,l,data(i,k),data(i,l));
    %                         result_theta(k,l,data(i,l),data(i,k)) = ...
    %                             result_theta(l,k,data(i,k),data(i,l));
                       end
                    end

                 end
        end
    end
end
theta1(f == 0)=1;
f = reshape(f,[numel(theta),1]);
theta1 = reshape(theta1,[numel(theta),1]);
%theta=linprog(f*-1,[f'; ones(1,numel(theta))*-1],[1 0],theta1',0,[],[]);
theta=linprog(f*-1, ones(1,numel(theta)), 0,theta1',0,ones(1,numel(theta))*-0,ones(1,numel(theta))*2);
if size(theta,1)==0
else
theta=reshape(theta,theta_shape);
end
end