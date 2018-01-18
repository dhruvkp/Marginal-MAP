% this is a simple BFS to create a tree out of Adjacent matrix
% where first column represent parents
% and secon represents childs
function input_tree=create_tree(adjac_mat)
input_tree=cell(size(adjac_mat,1),3);
visted=[1];
explored=[1];
while size(explored,2)>=1
    node_no=explored(1);
    [output_tree, visit, childs]=visit_node_childs(node_no,input_tree,adjac_mat,visted);
    explored=[explored childs];
    input_tree=output_tree;
    if size(explored,2)<1
    else
        explored=explored(:,2:size(explored,2));
    end
    visted=visit;
end
end