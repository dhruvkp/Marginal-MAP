function product = calc_msg_prod(input_tree)

% for each node calculate the mse product and store it in prod
product = cell(size(input_tree,1),1);

for i = 1:size(input_tree,1)
    if ~ isempty(input_tree{i, 1})
        if ~ isempty(product{i,1})
            product{i,1} = product{i,1} .* prod(input_tree{i, 1}(2:size(input_tree{i, 1},1)), 2);
        else
            product{i,1} = prod(input_tree{i, 1}(2:size(input_tree{i, 1},1)), 2);
        end
    end
    if ~ isempty(input_tree{i, 2})
        if ~ isempty(product{i,1})
            product{i,1} = product{i,1} .* prod(input_tree{i, 2}(2:size(input_tree{i, 2},1)), 2);
        else
            product{i,1} = prod(input_tree{i, 2}(2:size(input_tree{i, 2},1)), 2);
        end
    end
end
end