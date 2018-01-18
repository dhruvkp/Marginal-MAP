function [output_tree, visit, childs, out_theta]=visit_node_childs1(node_no,theta,input_tree,adjac_mat,visted)
childs=[];
visit=visted;
out_theta= theta;
for i=1:size(adjac_mat,1)
    if find(adjac_mat(node_no,i)==1 && size(find(visted==i),2)<1)
        visit=[visit i];
        input_tree{i,1}=[node_no;1 ;1 ];
        input_tree{node_no, 2}=[ input_tree{node_no,2} [i;1;1]];
        childs=[childs i];
        out_theta{node_no,i}=rand(input_tree{node_no,3},input_tree{i,3})*5+5;
    end
end
output_tree=input_tree;
end