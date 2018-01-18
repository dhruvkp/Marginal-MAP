function delta_theta=M_step1(trees,q_of_h,theta,indexes,training_data,par)


delta_theta = theta;
for i=1:size(training_data,1)
    temp_example=training_data{i,:};
    temp_example1=training_data{i,:};
    index = find([q_of_h{:,1}] == size(training_data{i,:},2));
    indexes.input = size(training_data{i,:},2);
    x = 1:size(training_data{i,:},2)-1;
    indexes.hidden = x(rem(x,2)==1);
    indexes.output = x(rem(x,2)==0);
    input_tree = trees{index , 2};
    temp_theta.h_to_y = zeros(par.values(1,1),par.values(2,1));
    temp_theta.h_to_h = zeros(par.values(1,1),par.values(1,1));
    temp_theta.h = zeros(par.values(1,1),1);
    for parent = 1:size(training_data{i,:},2)
        if isempty(input_tree{parent,2})
            continue;
        end
        for child = input_tree{parent,2}(1,:)
            temp_example1(1,indexes.hidden)=ones(size(indexes.hidden))*-1;
            temp_example1(1,indexes.output)=ones(size(indexes.output))*-1;
            p_of_observed=inference1(input_tree,theta,temp_example1,indexes);
            for parent_value = 1: input_tree{parent,3}
                for child_value = 1: input_tree{child,3}
                     
                        if ismember(child,indexes.input) && ismember(parent ,indexes.hidden)
                            temp_theta.h(parent_value ,child_value) = temp_theta.h(parent_value ,child_value) ...
                               +par.lRate ;
                        elseif ismember(parent,indexes.hidden) && ismember(child,indexes.hidden)
                            temp_theta.h_to_h(parent_value ,child_value) = temp_theta.h_to_h(parent_value ,child_value) ...
                                +par.lRate;
                        else
                            if temp_example(1,child ) == child_value
                                temp_theta.h_to_y(parent_value ,child_value) = temp_theta.h_to_y(parent_value ,child_value) ...
                                +par.lRate;
                            end
                        end
            
                               
                           
                    
                    
                    if (ismember(parent,indexes.input)&&temp_example1(1,parent)~=parent_value)...
                                        || (ismember(child,indexes.input)&& temp_example1(1,child)~=child_value)
                    else

                        temp_example1(1,parent)=parent_value;
                        temp_example1(1,child)=child_value;
                        dif_of_pobserved=inference1(input_tree,theta,temp_example1,indexes);
                        if ismember(parent,indexes.input) && ismember(child,indexes.hidden)
                            temp_theta.h(parent_value ,child_value) = temp_theta.h(parent_value ,child_value) ...
                               -(par.lRate*dif_of_pobserved/p_of_observed) ;
                        elseif ismember(parent,indexes.hidden) && ismember(child,indexes.hidden)
                            temp_theta.h_to_h(parent_value ,child_value) = temp_theta.h_to_h(parent_value ,child_value) ...
                                -(par.lRate*dif_of_pobserved/p_of_observed);
                        else
                            temp_theta.h_to_y(parent_value ,child_value) = temp_theta.h_to_y(parent_value ,child_value) ...
                                -(par.lRate*dif_of_pobserved/p_of_observed);
                        end
                        
                        
                        
                    end
                    
                end
            end
        end
    end
    
    delta_theta.h_to_y = delta_theta.h_to_y + temp_theta.h_to_y /size(temp_example,2) ;
    delta_theta.h_to_h = delta_theta.h_to_h + temp_theta.h_to_h / (size(temp_example,2) -1);
    delta_theta.h = delta_theta.h + temp_theta.h; 
    
    
end


delta_theta.h_to_y = delta_theta.h_to_y - 2* par.lambda * par.lRate *theta.h_to_y;
delta_theta.h_to_h = delta_theta.h_to_h - 2* par.lambda * par.lRate *theta.h_to_h;
delta_theta.h = delta_theta.h - 2* par.lambda * par.lRate *theta.h;


%delta_theta = delta_theta - 2 * par.lambda* theta;
end