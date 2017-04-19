A=[0 1 1;1 0 0;1 0 0]; %adjacency matrix

input_indexes=[1];
hidden_indexes=[2];
output_indexes=[3];

number_variables=size(A,1);
theta_c=zeros(3,3,2,2); %emptry parameters
x_domain=[1,2]; % binary variables domain

training_size=20;
test_size=20;

% gibbs sampling
burnin=10;
number_samples=220;
samples=gibbs_sampler_mrf_with_edge_parameters(A,theta_c,x_domain,burnin,training_size+test_size);
training_data=samples(1:training_size,:);
test_data=samples(1:test_size,:);
[theta_c] = EM( A,input_indexes,output_indexes,hidden_indexes,x_domain,training_data);