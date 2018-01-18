load('../trained_5000x5_words.mat');
addpath('../');
test_data_size=1000;
test_samples=zeros(test_data_size,2*max_line_size);
test_samples(:,output)=training_data(sample_size-test_data_size+(1:test_data_size),1:max_line_size);
test_samples(:,hidden)=-1;
fid = fopen('results.txt', 'w') ;
results={[],[],[],[],[]};
theta_i_n=ones(size(A,1),max(domain_sizes));
theta_i_n(1,:)=theta_i;
for i=1:100
    positions=2*randsample(max_line_size,characters_to_predict)'; % generate three random character positions to predict
    input=setdiff(output,positions);
    sample=test_samples(i,:);
    true_sample=sample;
    sample(positions)=-1;
    sample(positions)=predict_characters( A,positions,input,hidden,domain_sizes,theta_i_n,theta_c_learned,sample );
    fprintf(fid,'%d.\n%s\n%s\n\n',i,get_sentence(dictionary,true_sample,output),get_sentence(dictionary,sample,output));
    results{positions/2}=[results{positions/2},sample(positions)];
end
fclose(fid);