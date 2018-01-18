addpath('../../EM/V8')
% load('trained_5000x100.mat')
% load('../../datasets/MSR_sentence_completion/5000_words.mat');
% 
max_line_size=size(training_data,2);
hidden_domain_size=10;
algorithm=1; % 1 for FW, 2 for EM

%give row vectors
input=[];
output=2*(1:max_line_size);
hidden=2*(1:max_line_size)-1;
domain_sizes=zeros(2*max_line_size,1);
domain_sizes(output)=length(dictionary)*ones(length(output),1);
domain_sizes(hidden)=hidden_domain_size*ones(length(hidden),1);
A=zeros(2*max_line_size,2*max_line_size);
for i=1:max_line_size
    j=2*i-1;
    A(j,j+1)=1;
    A(j+1,j)=1;
    if j+2<2*max_line_size
        A(j,j+2)=1;
        A(j+2,j)=1;
    end
end

sample_size=4000;
samples=zeros(sample_size,2*max_line_size);
samples(:,output)=training_data(1:sample_size,1:max_line_size);
samples(:,hidden)=-1;

training_data_em=cell(size(samples,1),1);
for i=1:size(samples,1)
        training_data_em{i}=[samples(i,2:2:(2*max_line_size)),1];
end
% 
% theta_c_em=EM1(theta,param,indices,training_data);
theta_i=[];
theta_c_learned=[];
parameters.iterations=100;
parameters.tolerance=0.00001;
LR_range=[0.001,.01,.1,1,10,100];
% for LR_temp=LR
lambda_range=[0.01,0.1,1,10,100];
max_accuracy=0;
max_theta_i=[];
max_lr=0;
max_lambda=0;
max_theta_c=[];
for lr=LR_range
    for lambda=lambda_range
        parameters.LR=lr;
        parameters.lambda=lambda;
        
        if algorithm==1
            [theta_i,theta_c_learned]=learn(A,input,output,hidden,domain_sizes,samples,lambda,'linprog1');
        else
            theta_c_em=hmm_em(hidden_domain_size,length(dictionary),training_data_em,parameters);
            number_variables=size(A,1);
            theta_c_learned=zeros(number_variables,number_variables,max(domain_sizes),max(domain_sizes));
            theta_i=theta_c_em.S;
            for i=hidden
                theta_c_learned(i,i+1,1:hidden_domain_size,:)=theta_c_em.B;
                if i+2>number_variables
                    break
                end
                theta_c_learned(i,i+2,1:hidden_domain_size,1:hidden_domain_size)=theta_c_em.A;
            end
        end
        
        addpath('./inference/');
        addpath('./permn/');
        
        test_data_size=1000;
        test_samples=zeros(test_data_size,2*max_line_size);
        test_samples(:,output)=training_data(sample_size+(1:test_data_size),1:max_line_size);
        test_samples(:,hidden)=-1;
        fid = fopen('results.txt', 'w') ;
        characters_to_predict=2;
        total_correct=0;
        for i=1:test_data_size
            positions=2*randsample(max_line_size,characters_to_predict)'; % generate three random character positions to predict
            input=setdiff(output,positions);
            sample=test_samples(i,:);
            true_sample=sample;
            sample(positions)=-1;
            x=bp_predict_characters( A,output,[],hidden,positions,domain_sizes,theta_i,theta_c_learned,sample );
            sample(positions)=x;
            fprintf(fid,'%d.\n%s\n%s\n\n',i,get_sentence(dictionary,true_sample,output),get_sentence(dictionary,sample,output));
            total_correct=total_correct+characters_to_predict-nnz(true_sample-sample);
        end
        fclose(fid);
%         parameters.lambda
        accuracy=total_correct/(characters_to_predict*test_data_size)
        if accuracy>max_accuracy
            max_accuracy=accuracy;
            max_theta_i=theta_i;
            max_theta_c=theta_c_learned;
            max_lambda=lambda;
            max_lr=lr;
        end
    end
end