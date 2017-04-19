% calculate the probability of a particular example
function prob=calc_p_all_x(A,training_example,theta_c)
prob=1;
number_variables=size(A,1);
% for all edges in the graph multibly it's factor to the probability
for i=1:number_variables
    for j=i+1:number_variables
        if A(i,j)>0 
            prob=prob*exp(theta_c(i,j,training_example(1,i),training_example(1,j)));
       end
    end
end
end