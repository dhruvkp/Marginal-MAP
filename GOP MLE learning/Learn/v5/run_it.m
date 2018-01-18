addpath('../QPBO-v1.4.src')
addpath('../../EM/V5')


% A=[0 1 1;1 0 0;1 0 0]; %adjacency matrix
% hidden=[1]; %#ok<*NBRAK>
% input=[2];
% output=[3];
% sample_sizes=[200,400];


A= [0 1 1 0 0 0;1 0 0 0 0 0;1 0 0 1 1 0;0 0 1 0 0 0;0 0 1 0 0 1;0 0 0 0 1 0];

%give row vectors
input=[2,4];
hidden=[1,3,5];
output=[6];
sample_sizes=[100,200,400,600,800];
domain_sizes=[6,7,7,2,3,5];

% A=[0 1;1 0];
% 
% %give row vectors
% input=[1];
% hidden=[];
% output=[2];
% sample_sizes=[40,60,80];

number_runs=10;


number_variables=size(A,1);
max_domain_size=max(domain_sizes);

theta_c=rand(number_variables,number_variables,max_domain_size,max_domain_size)*10-5;
for i=1:number_variables
    for j=i:number_variables
        if A(i,j)==0
            theta_c(i,j,:,:)=zeros(max_domain_size,max_domain_size);
            theta_c(j,i,:,:)=theta_c(i,j,:,:);
        end
        theta_c(i,j,:,:)=squeeze(theta_c(j,i,:,:))';
    end
end

KLD_values_qpbo=zeros(number_runs,length(sample_sizes));
KLD_values=zeros(number_runs,length(sample_sizes));
KLD_values_em=zeros(number_runs,length(sample_sizes));
status=zeros(1,length(sample_sizes));
cnt=1;
burnin=100;

for training_size=sample_sizes
    % gibbs sampling
    %
    % % initial accuracy
    % theta_c
    % [ out_list2, out_prob2 ]=generate_output( A,input,output,hidden,x_domain,test_data,theta_c);
    % acc_before=calc_accuracy(out_list2,test_data,output)
    samples=gibbs_sampler_mrf_with_edge_parameters(A,theta_c,domain_sizes,burnin,training_size);
    
    
%     for n=1:number_runs
%         lambda=lambda_range(n);
%         samples=gibbs_sampler_mrf_with_edge_parameters(A,theta_c,x_domain,burnin,training_size);
        training_data=samples(1:training_size,:);
        %test_data=samples(1:test_size,:);
        
        % remove hidden variables from training_data
        training_data(:,hidden)=0;
        %learn parameters using structured learning
        %tic
        %learn_merged(A,input,output,hidden,x_domain,training_data,1);
        %toc
%         tic

%          theta_c_learned=learn(A,input,output,hidden,x_domain,training_data,1,'compare');

        theta_c_learned_qpbo=learn(A,input,output,hidden,[1,2],training_data,1,'qpbo');
%        theta_c_learned_em=learn(A,input,output,hidden,x_domain,training_data,1,'linprog2');
%        theta_c_learned=learn(A,input,output,hidden,domain_sizes,training_data,1,'linprog1');


         

%         toc
%         tic
        [theta_c_em] = EM1( A,domain_sizes,hidden,output,input,training_data,2);
%         toc
%          theta_c_learned_all(n,:,:,:,:)=theta_c_learned;
%          theta_c_learned_all_qpbo(n,:,:,:,:)=theta_c_learned_qpbo;
        
%         KLD_values(n,cnt)=KL_divergence( A,output,input,hidden,theta_c,theta_c_learned);
%         KLD_values_qpbo(n,cnt)=KL_divergence( A,output,input,hidden,theta_c,theta_c_learned_qpbo);
%         KLD_values_em(n,cnt)=KL_divergence( A,output,input,hidden,theta_c,theta_c_em);
%         [~,KLD_values_qpbo(:,cnt)]=parameter_tuning( A,x_domain,input,output,hidden,training_data,theta_c,number_runs,'qpbo' );
%         [~,KLD_values_em(:,cnt)]=parameter_tuning( A,x_domain,input,output,hidden,training_data,theta_c,number_runs,'EM' );
        [cnt;training_size]
    
%     end
    cnt=cnt+1;
end
%
figure
errorbar(sample_sizes,mean(KLD_values,1),min(KLD_values,[],1)-mean(KLD_values,1),max(KLD_values,[],1)-mean(KLD_values,1),'-g.')
% hold on
% errorbar(sample_sizes,mean(KLD_values_qpbo,1),min(KLD_values_qpbo,[],1)-mean(KLD_values_qpbo,1),max(KLD_values_qpbo,[],1)-mean(KLD_values_qpbo,1),'-b.')
% hold on
% errorbar(sample_sizes,mean(KLD_values_em,1),min(KLD_values_em,[],1)-mean(KLD_values_em,1),max(KLD_values_em,[],1)-mean(KLD_values_em,1),'-r.')
%plot(sample_sizes,KLD_values,'-og',sample_sizes,KLD_values_em,'-ok');
ylabel('KL divergence of p(Y|X,theta)');
xlabel('Sample sizes');
legend('qpbo','EM')