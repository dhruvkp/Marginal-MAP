function q_of_h=E_step1(q_of_h,input_tree,theta,hidden_indexes,training_data)

for i=1:size(training_data,1)
    temp_example=training_data(i,:);
    
    for j=1:size(q_of_h,1)
        
        % add current h to the training example
        temp_example(1,hidden_indexes)=q_of_h(j,1:size(hidden_indexes,2));
        % calc P(h | x,y)
        q_of_h(j,i+size(hidden_indexes,2))= ...
            (inference1(input_tree,theta,temp_example)/...
                  inference1(input_tree,theta,training_data(i,:)));
    end
end
end