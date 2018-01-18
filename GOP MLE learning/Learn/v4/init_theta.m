function [ theta_i,theta_c ] = init_theta( A,domain_size )
%INIT_THETA Summary of this function goes here
%   Detailed explanation goes here
number_variables=size(A,1);
number_edges=sum(sum(A))/2;

theta_i=10*rand(number_variables,domain_size);
theta_c=10*rand(number_edges,domain_size,domain_size);

end

