function [ z ] = compute_exact_z( A,theta_c )
%COMPUTE_EXACT_Z computes partition function
number_variables=size(A,1);
all_assignments=binary_combinations(number_variables);
z=0;
for i=all_assignments'
    x=i';
    z=z+calc_p_all_x(A,x,theta_c);
end
end

