function tree1=add_childs(adjac_mat, parent_id, tree)

disp(tree.tostring)
    for j=1:size(adjac_mat,1)
        if(adjac_mat(parent_id,j)==1)
            [ tree1 node2 ] = tree.addnode(parent_id, num2str(1));
            adjac_mat(j,parent_id)=0;
        end

    end
end