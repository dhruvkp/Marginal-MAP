function [] = dataset_init()
%MSR_INIT Summary of this function goes here
%   Detailed explanation goes here

training_data_dir='Training_Data';
% hidden_variables_domain_size=36; %number of possible POS tags are 36 as per https://www.ling.upenn.edu/courses/Fall_2003/ling001/penn_treebank_pos.html
fnames=dir(training_data_dir);
fnames=fnames(3:end);
dictionary={};
max_line_size=50;
max_word_size=5;
samples=5000;
training_data=zeros(samples,max_word_size);
line_number=0;
word_number=0;
for cnt=1:length(fnames)
    f=fnames(cnt);
    f.name
    text=fileread(strcat(training_data_dir,'\',f.name));
    text(text==char(10))=' ';
    text=regexp(text,'(?<=[.!?])\s* ','split');
    lengths=[];
    for line=text(75:end)
        words=strsplit(line{1},'[^a-zA-Z\d'']','DelimiterType','RegularExpression');
        for word=words
            if length(word{1})~=max_word_size
                continue
            end
            lower_word=lower(word{1});
            sample=[];
            for k=1:max_word_size
                [t,i]=ismember(lower_word(k),dictionary);
                if ~t
                    dictionary=[dictionary,lower_word(k)];
                    i=length(dictionary);
                end
                sample=[sample,i];
            end
            if sum(ismember(training_data, sample, 'rows'))
                continue
            end
            word_number=word_number+1;
            training_data(word_number,:)=sample;
            if word_number==samples
                break
            end
        end
        if word_number==samples
            break
        end
%         line_array=line{1};
%         k=0;
%         if length(line_array)<max_line_size
%             continue
%         end
%         line_number=line_number+1;
% 
%         for k=1:length(line_array)
%             if k>max_line_size
%                 break
%             end
%             lower_char=lower(line_array(k));
%             [t,i]=ismember(lower_char,dictionary);
%             if ~t
%                 dictionary=[dictionary,lower_char];
%                 i=length(dictionary);
%             end
%             if i==0
%             end
%             training_data(line_number,k)=i;
%         end
%         if k<max_line_size
%             training_data(line_number,(k+1):max_line_size)=ones(1,max_line_size-k);
%         end
%         if line_number>=samples
%             break
%         end
    end
    if word_number==samples
        break
    end
    if line_number>=samples
        break
    end
end
save(strcat(num2str(samples),'_words.mat'),'training_data','dictionary');
end

