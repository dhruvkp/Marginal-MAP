load('../../datasets/MSR_sentence_completion/5000_words.mat');
max_line_size=size(training_data,2);
hidden_domain_size=10;

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
sample_size=5000;
samples=zeros(sample_size,2*max_line_size);
samples(:,output)=training_data(1:sample_size,1:max_line_size);
samples(:,hidden)=-1;
%     [theta_i,theta_c_learned]=learn(A,input,output,hidden,domain_sizes,samples,1,'linprog1');
% 
for lambda=[0.001,0.01,0.1]
    [theta_i,theta_c_learned]=learn(A,input,output,hidden,domain_sizes,samples,lambda,'linprog1');
    save(strcat('trained_5000x5_words_',string(lambda)));
end