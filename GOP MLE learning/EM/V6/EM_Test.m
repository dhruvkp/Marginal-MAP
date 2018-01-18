% adjac_mat=[0 1 1;1 0 0;1 0 0]; %adjacency matrix
% 
% indexes.input=[2];
% indexes.hidden=[1];
% indexes.output=[3];
adjac_mat=[0 1 1 0 0 0 1;
    1 0 0 0 0 0 0;
    1 0 0 1 1 0 0;
    0 0 1 0 0 0 0;
    0 0 1 0 0 1 0;
    0 0 0 0 1 0 0;
    1 0 0 0 0 0 0];

%give row vectors
indexes.input=[7];
indexes.hidden=[1,3,5];
indexes.output=[2,4,6];
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
training_data(:,indexes.hidden)=-1;
test_data=samples(1:test_size,:);
training_data(:,indexes.hidden)=-1;
training_data(:,indexes.input)=-1;

clear theta;
%theta=EM(adjac_mat,ones(number_variables,1)*2,size(theta),indexes.hidden,indexes.output,indexes.input,training_data,2);
values = ones(2,1)*4;
param.values=values;
param.lambda=2;
param.lRate = .01;
param.threshold = .001;
param.iterations = 50;
%input_tree =create_tree1(adjac_mat , param);
theta.h_to_y = rand(4,4);
theta.h_to_h = rand(4,4);
theta.h = rand(4,1)*5+5;
clear training_data;
training_data{1,1}=[-1 2 1];
training_data{2,1}=[-1 2 -1 4 1];
theta1=EM1(theta,param,indexes,training_data);