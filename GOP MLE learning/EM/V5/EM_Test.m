% adjac_mat=[0 1 1;1 0 0;1 0 0]; %adjacency matrix
% 
% input_indexes=[2];
% hidden_indexes=[1];
% output_indexes=[3];
adjac_mat=[0 1 1 0 0 0;
    1 0 0 0 0 0;
    1 0 0 1 1 0;
    0 0 1 0 0 0;
    0 0 1 0 0 1;
    0 0 0 0 1 0];

%give row vectors
input_indexes=[2,4];
hidden_indexes=[1,3,5];
output_indexes=[6];
x_domain=[1,2]; % binary variables domain
number_variables=size(adjac_mat,1);
theta=zeros(number_variables,number_variables,2,2); %emptry parameters
for i=1:number_variables
    for j=1:number_variables
        if adjac_mat(i,j)>0
            for x_i=x_domain
                for x_j=x_domain
                    theta(i,j,x_i,x_j)=(x_i+x_j);
                end
            end
        end
    end
end


training_size=10;
test_size=10;

% gibbs sampling
burnin=1000;
number_samples=2;
samples=gibbs_sampler_mrf_with_edge_parameters(adjac_mat,theta,x_domain,burnin,training_size+test_size);
training_data=samples(1:training_size,:);
training_data(:,hidden_indexes)=-1;
test_data=samples(1:test_size,:);
%theta=EM(adjac_mat,ones(number_variables,1)*2,size(theta),hidden_indexes,output_indexes,input_indexes,training_data,2);
theta1=EM1(adjac_mat,ones(number_variables,1)*4,hidden_indexes,output_indexes,input_indexes,training_data,2);