function theta=EM1(theta , par,indexes,training_data)
% par.values is a list of maximum value any node can get
% mark all hidden values as a blank
% Create a tree out of our graph
[q_of_h , trees ] = create_structures(training_data , par.values)

temp_theta=theta.h_to_y+1;
count=0;
threshold=ones(size(theta.h_to_y))* par.threshold;
result1=abs(temp_theta-theta.h_to_y)>threshold;
while  size(find(result1==1),1)>0 && count < par.iterations
    result1=abs(temp_theta-theta.h_to_y)>threshold;
    temp_theta=theta.h_to_y;
    count=count+1;
   % q_of_h=E_step1(q_of_h,input_tree,theta,indexes,training_data);
    theta=M_step1(trees,q_of_h,theta,indexes,training_data,par);
end
end

function [q_of_h , trees ] = create_structures(training_data , values)
% get all different sizes of the training data
cellsz = unique([cellfun(@length,training_data)],'rows');
% q_of_h is a cell of all combination of h's, that will be required to
% train out model
trees = cell(size(cellsz,1),1);
q_of_h = cell(size(cellsz,1),1);
q_of_h= [ mat2cell(cellsz,ones(size(cellsz)),1) , q_of_h ];
trees= [ mat2cell(cellsz,ones(size(cellsz)),1) , trees ];
% create a table out of all posible combinations of h and add
%columns for the reuslts of every training example

for i=1:size(cellsz,1)
    q_of_h{i,2}=combinator(values(1,1),(q_of_h{i,1}-1)/2,'p','r');
    all_values = ones(trees{i,1},1);
    all_values(2:2:length(all_values)) = values(2,1) ;
    all_values(1:2:length(all_values)) = values(1,1) ;
    all_values(trees{i,1},1)=1;
    trees{i,2}=create_tree1(trees{i,1}, all_values);
    q_of_h{i,2}=[ q_of_h{i,2}, zeros(size(q_of_h{i,2},1),1)];
end
end