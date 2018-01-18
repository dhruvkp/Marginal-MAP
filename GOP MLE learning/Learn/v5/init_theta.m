function [ theta_i,theta_c ] = init_theta( A,domain_sizes )
%INIT_THETA Summary of this function goes here
%   Detailed explanation goes here
number_variables=size(A,1);
number_edges=sum(sum(A))/2;
domain_size=max(domain_sizes);

theta_i=10*ones(number_variables,domain_size);
theta_c=10*ones(number_edges,domain_size,domain_size);

end

