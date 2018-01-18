% Input: given theta and the graph information 
% output: Expected out put values
function [out_list , out_prob]=generate_output( A,input_indexes,output_indexes,hidden_indexes,x_domain,training_data,theta_c)
number_variables=size(A,1);
data_size=size(training_data,1);
% all possible values for hte hidden variables combined
hidden_vals=binary_combinations(size(hidden_indexes,1));
% all possible values for the output values combined
out_vals=binary_combinations(size(output_indexes,1));
out_list=zeros(data_size,size(output_indexes,1));
out_prob=zeros(data_size,size(out_vals,1));


for i=1:data_size
    max=-inf;
    for j=1:size(out_vals,1)
        prob=0;
        training_data(i,output_indexes)=out_vals(j,:);
        for k=1:size(hidden_vals,1)
            training_data(i,hidden_indexes)=hidden_vals(k,:);
            prob=prob+calc_p_all_x(A,training_data(i,:),theta_c);
        end
        out_prob(i,j)=prob;
        if prob >= max
            out_list(i,:)=out_vals(j,:);
            max=prob;
        end
    end
    
end

end