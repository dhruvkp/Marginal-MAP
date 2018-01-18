% A=[0 1 1;1 0 0;1 0 0]; %adjacency matrix
% hidden=[];
% input=[1,2];
% output=[3];
A=[0 1 1 0 0 0;1 0 0 0 0 0;1 0 0 1 1 0;0 0 1 0 0 0;0 0 1 0 0 1;0 0 0 0 1 0];
%A=[0 1;1 0]
% give row vectors
input=[2,4];
hidden=[1,3,5];
output=[6];

number_variables=size(A,1);
theta_c=zeros(number_variables,number_variables,2,2); %emptry parameters
x_domain=[1,2]; % binary variables domain
domain_size=length(x_domain);

for i=1:number_variables
    for j=1:number_variables
        if A(i,j)>0
            for x_i=x_domain
                for x_j=x_domain
                    theta_c(i,j,x_i,x_j)=(x_i+x_j)*10;
                end
            end
        end
    end
end

number_runs=1;
sample_sizes=[700];
theta_c_learned_all=zeros(number_runs,number_variables,number_variables,domain_size,domain_size);
KLD_values=zeros(number_runs,length(sample_sizes));
KLD_values_em=zeros(number_runs,length(sample_sizes));
status=zeros(1,length(sample_sizes));
cnt=1;
for training_size=sample_sizes
    % gibbs sampling
    burnin=100;
    %
    % % initial accuracy
    % theta_c
    % [ out_list2, out_prob2 ]=generate_output( A,input,output,hidden,x_domain,test_data,theta_c);
    % acc_before=calc_accuracy(out_list2,test_data,output)
    
    
    for n=1:number_runs
        
        samples=gibbs_sampler_mrf_with_edge_parameters(A,theta_c,x_domain,burnin,training_size);
        training_data=samples(1:training_size,:);
        %test_data=samples(1:test_size,:);
        
        % remove hidden variables from training_data
        training_data(:,hidden)=0;
        %learn parameters using structured learning
        %theta_c_learned=learn(A,input,output,hidden,x_domain,training_data,1);
        [theta_c_em] = EM( A,input,output,hidden,x_domain,training_data);
        theta_c_learned_all(n,:,:,:,:)=theta_c_learned;
        
        %KLD_values(n,cnt)=KL_divergence( A,output,input,hidden,theta_c,theta_c_learned);
        KLD_values_em(n,cnt)=KL_divergence( A,output,input,hidden,theta_c,theta_c_em);

    
    end
    cnt=cnt+1;
end
%
figure
errorbar(sample_sizes,mean(KLD_values),min(KLD_values)-mean(KLD_values),max(KLD_values)-mean(KLD_values),'-g.')
hold on
errorbar(sample_sizes,mean(KLD_values_em),min(KLD_values_em)-mean(KLD_values_em),max(KLD_values_em)-mean(KLD_values_em),'-r.')
%plot(sample_sizes,KLD_values,'-og',sample_sizes,KLD_values_em,'-ok');
ylabel('KL divergence of p(Y|X,theta)');
xlabel('Number of samples');
legend('Frank-wolfe','EM algorithm')