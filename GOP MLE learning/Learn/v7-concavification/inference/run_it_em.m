load('../trained_5000x5_words_10h_EM.mat');
addpath('../');
theta_i=theta_c_em.S;
number_variables=size(A,1);
theta_c_learned=zeros(number_variables,number_variables,max(domain_sizes),max(domain_sizes));
for i=hidden
    if i+2>number_variables
        break
    end
    theta_c_learned(i,i+2,1:hidden_domain_size,1:hidden_domain_size)=theta_c_em.A;
    theta_c_learned(i,i+1,1:hidden_domain_size,:)=theta_c_em.B;
end
test_data_size=1000;
test_samples=samples(sample_size-test_data_size+(1:test_data_size),:);
fid = fopen('results.txt', 'w') ;
% results={[],[],[],[],[]};
characters_to_predict=2;
total_correct=0;
for i=1:test_data_size
positions=2*randsample(max_line_size,characters_to_predict)'; % generate three random character positions to predict
input=setdiff(output,positions);
sample=test_samples(i,:);
true_sample=sample;
sample(positions)=-1;
sample(positions)=bp_predict_characters( A,output,[],hidden,positions,domain_sizes,theta_i,theta_c_learned,sample );
fprintf(fid,'%d.\n%s\n%s\n\n',i,get_sentence(dictionary,true_sample,output),get_sentence(dictionary,sample,output));
% results{positions/2}=[results{positions/2},sample(positions)];
total_correct=total_correct+characters_to_predict-nnz(true_sample-sample);
end
fclose(fid);
accuracy=total_correct/(characters_to_predict*test_data_size)