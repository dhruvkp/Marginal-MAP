function theta=EM1(adjac_mat, node_value,hidden_indexes,output_indexes,input_indexes,training_data,lambda)

% mark all hidden values as a blank
training_data(:,hidden_indexes)=-1;
% Create a tree out of our graph
[ input_tree , theta ]=create_tree1(adjac_mat , node_value);



% create a table out of all posible combinations of h and add
%columns for the reuslts of every training example
q_of_h=combinator(input_tree{hidden_indexes(1,1),3},size(hidden_indexes,2),'p','r');
q_of_h=[ q_of_h, zeros(size(q_of_h,1),size(training_data,1))];

    q_of_h=E_step1(q_of_h,input_tree,theta,hidden_indexes,training_data);
    temp_q_of_h=q_of_h+1;
    count=0;
    threshold=ones(size(q_of_h))*.001;
    result1=abs(temp_q_of_h-q_of_h)>threshold;
    learn_rate=0.01;
while  size(find(result1==1),1)>0 && count<50
    result1=abs(temp_q_of_h-q_of_h)>threshold;
    temp_q_of_h=q_of_h;
    count=count+1;
    theta=M_step1(adjac_mat,q_of_h,input_tree,theta,hidden_indexes,output_indexes,input_indexes,training_data,learn_rate,lambda);
    %theta=theta-2*lambda*theta;
   % theta=theta+0.01*delta_theta;
    q_of_h=E_step1(q_of_h,input_tree,theta,hidden_indexes,training_data);
end
end
