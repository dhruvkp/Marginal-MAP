A=[0 0 0 0 0 1 0;
   0 0 0 0 0 1 0;
   0 0 0 0 0 1 0;
   0 0 0 0 0 1 0;
   0 0 0 0 0 1 0;
   1 1 1 1 1 0 1;
   0 0 0 0 0 1 0;]; %adjacency matrix

input_indexes=[1 ;2; 3; 4; 5];
hidden_indexes=[6];
output_indexes=[7];


number_variables=size(A,1);
% binary variables domain
x_domain=[1,2];
% initialize_theta
theta_c=initialize_theta(A,number_variables,x_domain); %emptry parameters
 
% select data sizes
training_size=250;
test_size=50;

% gibbs sampling
burnin=10000;
number_samples=20;
samples=gibbs_sampler_mrf_with_edge_parameters(A,theta_c,x_domain,burnin,training_size+test_size);
training_data=samples(1:training_size,:);
test_data=samples(training_size:test_size,:);

% measuring accuracy using the original theta
[ out_list1, out_prob1]=generate_output( A,input_indexes,output_indexes,hidden_indexes,x_domain,training_data,theta_c);
acc_befor=calc_accuracy(out_list1,training_data,output_indexes);

% train the modle
[theta_c] = EM( A,input_indexes,output_indexes,hidden_indexes,x_domain,training_data);
[ out_list2, out_prob2 ]=generate_output( A,input_indexes,output_indexes,hidden_indexes,x_domain,training_data,theta_c);

% measuring accuracy using the trained theta
acc_after=calc_accuracy(out_list2,training_data,output_indexes);

% calculating the ratio between p(y,x) before and after training
temp_p1=sum(out_prob1,2);
temp_p2=sum(out_prob2,2);
temp_p1=repmat(temp_p1,1,2);
temp_p2=repmat(temp_p2,1,2);
out_prob1=out_prob1 ./ temp_p1;
out_prob2=out_prob2 ./ temp_p2;

probs=out_prob1-out_prob2;
