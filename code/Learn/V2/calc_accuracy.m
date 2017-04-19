function accuracy=calc_accuracy(out_list,training_data,output_indexes)
accuracy=0;
% count no. of orginal outputs that match the expected output
for i=1:size(out_list,1)
    if out_list(i,:)== training_data(i,output_indexes)
        accuracy=accuracy+1;
    end
end
accuracy=accuracy/size(out_list,1);
end