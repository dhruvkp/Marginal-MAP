% this is a simple BFS to create a tree out of Adjacent matrix
% where first column represent parents
% and second represents childs
function [ input_tree , theta, lambda ] =create_tree(adjac_mat , node_value )
theta = cell(size(adjac_mat,1),size(adjac_mat,1));
input_tree=[cell(size(adjac_mat,1),2), mat2cell(node_value,ones(size(node_value)),1)];
visted=[1];
explored=[1];
%input_tree{:,3}= ;
while size(explored,2)>=1
    node_no=explored(1);
    
    [output_tree, visit, childs, theta]=visit_node_childs(node_no,theta,input_tree,adjac_mat,visted);
    explored=[explored childs];
    input_tree=output_tree;
    if size(explored,2)<1
    else
        explored=explored(:,2:size(explored,2));
    end
    visted=visit;
end
lambda = max(sum(adjac_mat,2))*2  - 1 ;
end