% calculate the probability of a particular example
function prob=calc_p_all_x(A,training_example,theta_i,theta_c)
prob=1;
number_variables=size(A,1);
% for all edges in the graph multibly it's factor to the probability
for i=1:number_variables
    prob=prob*exp(theta_i(i,training_example(i)));
    for j=i+1:number_variables
        if A(i,j)>0 
            prob=prob*exp(theta_c(i,j,training_example(i),training_example(j)));
       end
    end
end
end