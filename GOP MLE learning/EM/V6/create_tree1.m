% this is a simple BFS to create a tree out of Adjacent matrix
% where first column represent parents
% and secon represents childs
function input_tree  =create_tree1(variable_no , values )

visted=[1];
explored=[1];
adjac_mat = calc_mat(variable_no);
input_tree=[cell(size(adjac_mat,1),2), mat2cell(values,ones(size(values)),1)];
%input_tree{:,3}= ;
while size(explored,2)>=1
    node_no=explored(1);
    
    [output_tree, visit, childs]=visit_node_childs1(node_no,input_tree,adjac_mat,visted);
    explored=[explored childs];
    input_tree=output_tree;
    if size(explored,2)<1
    else
        explored=explored(:,2:size(explored,2));
    end
    visted=visit;
    
end
end
function adjac_mat = calc_mat(variable_no)
adjac_mat = zeros(variable_no , variable_no);
for i = 1: variable_no-1
    if rem(variable_no,i+2)==0
        adjac_mat(variable_no , 1)=1;
        adjac_mat(1 , variable_no )=1;
    elseif rem(i,2)==0
        adjac_mat(i , i-1)=1;
        adjac_mat(i-1 , i )=1;
    else
        adjac_mat(i , i+2)=1;
        adjac_mat(i+2 , i )=1;
    end
end
end