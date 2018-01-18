function q_of_h=E_step(adjac_mat,input_tree,theta,hidden_indexes,training_data)
% create a table out of all posible combinations of h and add
%columns for the reuslts of every training example
q_of_h=binary_combinations(size(hidden_indexes,2));
q_of_h=[ q_of_h, zeros(size(q_of_h,1),size(training_data,1))];
for i=1:size(training_data,1)
    temp_example=training_data(i,:);
    
    for j=1:size(q_of_h,1)
        count=1;
        % add current h to the training example
        for k=1:size(hidden_indexes,2)
            temp_example(1,hidden_indexes(1,k))=q_of_h(j,count);
            count=count+1;
        end
        % calc P(h | x,y)
        q_of_h(j,i+size(hidden_indexes,2))= ...
            (calc_p_all_x(adjac_mat,temp_example,theta)/...
                  inference(input_tree,theta,training_data(i,:)));
    end
end
end