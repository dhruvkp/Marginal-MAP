function list=binary_combinations(att_no)
% att_no represent the number of attributes
list=de2bi(1:(2^att_no));
list=list(:,1:size(list,2)-1);
list(list==1)=2;
list(list==0)=1;
end