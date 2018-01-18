A=[0 1 1;1 0 0;1 0 0]; %adjacency matrix

% give row vectors
input=[2];
hidden=[1];
output=[3];

number_variables=size(A,1);
theta_c=zeros(3,3,2,2); %emptry parameters
x_domain=[1,2]; % binary variables domain
domain_size=length(x_domain);

for i=1:number_variables
    for j=1:number_variables
        if A(i,j)>0
            for x_i=x_domain
                for x_j=x_domain
                    theta_c(i,j,x_i,x_j)=(x_i+x_j)/12;
                end
            end
        end
    end
end

number_runs=2;
sample_sizes=[50];
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
    samples=gibbs_sampler_mrf_with_edge_parameters(A,theta_c,x_domain,burnin,training_size);
    training_data=samples(1:training_size,:);
    %test_data=samples(1:test_size,:);
    
    % remove hidden variables from training_data
    training_data(:,hidden)=0;
    
    for n=1:number_runs
        
        
        %learn parameters using structured learning
        theta_c_learned=learn(A,input,output,hidden,x_domain,training_data,50000);
        [theta_c_em] = EM( A,input,output,hidden,x_domain,training_data);
        theta_c_learned_all(n,:,:,:,:)=theta_c_learned;
        
        %%% KLD between p(Y|X,theta) and q(Y|X,theta) %%%
        
        % computing z %
        z1=0;
        z2=0;
        z3=0;
        for xi=x_domain
            for xj=x_domain
                prob1=0;
                prob2=0;
                prob3=0;
                for xk=x_domain
                    prob1=prob1+calc_p_all_x(A,[xk,xi,xj],theta_c);
                    prob2=prob2+calc_p_all_x(A,[xk,xi,xj],theta_c_learned);
                    prob3=prob3+calc_p_all_x(A,[xk,xi,xj],theta_c_em);
                end
                z1=z1+prob1;
                z2=z2+prob2;
                z3=z3+prob3;
            end
        end
        
        
        KL_divergence=0;
        KL_divergence_em=0;
        for xi=x_domain
            for xj=x_domain
                prob1=0;
                prob2=0;
                prob3=0;
                for xk=x_domain
                    prob1=prob1+calc_p_all_x(A,[xk,xi,xj],theta_c)/z1;
                    prob2=prob2+calc_p_all_x(A,[xk,xi,xj],theta_c_learned)/z2;
                    prob3=prob3+calc_p_all_x(A,[xk,xi,xj],theta_c_em)/z3;
                end
                KL_divergence=KL_divergence+prob1*log(prob1/prob2);
                KL_divergence_em=KL_divergence_em+prob1*log(prob1/prob3);
            end
        end
        KLD_values(n,cnt)=KL_divergence;
        KLD_values_em(n,cnt)=KL_divergence_em;
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
