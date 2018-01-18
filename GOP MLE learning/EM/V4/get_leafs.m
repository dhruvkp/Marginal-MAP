function leaf_list= get_leafs(input_tree)
leaf_list=[];
for i=1:size(input_tree,1)
    if size(input_tree{i,2},2)<1
        leaf_list=[leaf_list i];
    end
end
end