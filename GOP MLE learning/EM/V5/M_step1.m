function delta_theta=M_step1(adjac_mat,q_of_h,input_tree,theta,hidden_indexes,output_indexes,input_indexes,training_data,learn_rate,lambda)


delta_theta = theta;
for i=1:size(training_data,1)
    temp_example=training_data(i,:);
    temp_example1=training_data(i,:);
    
    for parent = 1:size(adjac_mat,1)
        if isempty(input_tree{parent,2})
            continue;
        end
        for child = input_tree{parent,2}(1,:)
            temp_example1(1,hidden_indexes)=ones(size(hidden_indexes))*-1;
            temp_example1(1,output_indexes)=ones(size(output_indexes))*-1;
            p_of_observed=inference1(input_tree,theta,temp_example1);
            for parent_value = 1: input_tree{parent,3}
                for child_value = 1: input_tree{child,3}
                    for hidden=1:size(q_of_h,1)
                        temp_example(1,hidden_indexes)=q_of_h(hidden,1);
                    
                            if temp_example(1,parent)== parent_value && temp_example(1,child)== child_value
                                delta_theta{parent , child}(parent_value ,child_value )= ...
                                    delta_theta{parent , child}(parent_value ,child_value )+learn_rate*q_of_h(hidden,i+size(hidden_indexes,2));
                               
                            end
                    end
                    
                    if (ismember(parent,input_indexes)&&temp_example1(1,parent)~=parent_value)...
                                        || (ismember(child,input_indexes)&& temp_example1(1,child)~=child_value)
                    else

                        temp_example1(1,parent)=parent_value;
                        temp_example1(1,child)=child_value;
                        dif_of_pobserved=inference1(input_tree,theta,temp_example1);
                        delta_theta{parent , child}(parent_value ,child_value )= ...
                                    delta_theta{parent , child}(parent_value ,child_value )-(learn_rate*dif_of_pobserved/p_of_observed);
                        
                    end
                    
                end
            end
        end
    end
    
    
    
    
    
   
end

for parent = 1:size(adjac_mat,1)
        if isempty(input_tree{parent,2})
            continue;
        end
        for child = input_tree{parent,2}(1,:)
            delta_theta{parent , child} = ...
                delta_theta{parent , child} - 2* lambda * learn_rate *theta{parent , child};
        end
end

%delta_theta = delta_theta - 2 * lambda* theta;
end