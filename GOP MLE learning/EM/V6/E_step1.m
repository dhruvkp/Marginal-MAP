function q_of_h=E_step1(q_of_h,input_tree,theta,indexes,training_example)

    temp_example=training_example;
    index = find([q_of_h{:,1}] == size(temp_example,2));
    for j=1:size(q_of_h{index , 2},1)
        
        % add current h to the training example
        temp_example(1,indexes.hidden)=q_of_h{index , 2}(j,1:size(indexes.hidden,2));
        % calc P(h | x,y)
        q_of_h{index , 2}(j,1+size(indexes.hidden,2))= ...
            (inference1(input_tree,theta,temp_example,indexes)/...
            inference1(input_tree,theta,training_example,indexes));
    end
end
