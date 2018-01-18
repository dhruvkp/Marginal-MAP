function [ predictions ] = predict_characters( A,output,input,hidden,domain_sizes,theta_i,theta_c,sample )
%PREDICT Summary of this function goes here
%   Detailed explanation goes here
% define order of variables to sum out
number_variables=size(A,1);
order_variables=horzcat(hidden,output,input);
domain_size=max(domain_sizes);

number_edges=sum(sum(A))/2;

% storing cliques
cliques=zeros(number_edges,2);
cnt=1;
for i=1:number_variables
    for j=i+1:number_variables
        if A(i,j)>0
            cliques(cnt,:)=[i j];
            cnt=cnt+1;
        end
    end
end


b_i=zeros(number_variables,domain_size);
b_c=zeros(number_edges,domain_size,domain_size);
for i=input
    b_i(i,sample(i))=1;
end
for i=setdiff(1:number_variables,input)
    b_i(i,:)=1/domain_sizes(i);
end
for c=1:number_edges
    i=cliques(c,1);
    j=cliques(c,2);
    if ismember(i,input) && ismember(j,input)
        b_c(c,sample(i),sample(j))=1;
    elseif ismember(i,input) && ~ismember(j,input)
        b_c(c,sample(i),:)=1/domain_sizes(j);
    elseif ~ismember(i,input) && ismember(j,input)
        b_c(c,:,sample(j))=1/domain_sizes(i);
    else
        b_c(c,:,:)=1/(domain_sizes(i)*domain_sizes(j));
    end
end

[ Aeq1,Aeq2,Beq1,Beq2,lb,ub ] = linprog_init( A,domain_sizes,cliques );

[w_i,w_i_c] = init_w( A,hidden );

%%%% learning with frank-wolfe %%%%

s_i=b_i;
s_c=b_c;

t=0;
while true
    t=t+1;
    [ dl_i,dl_c ] = calc_dl( b_i,b_c,theta_i,theta_c,w_i,w_i_c,cliques,order_variables,hidden,output,domain_sizes );
    gap=abs(sum(sum((b_i-s_i).*dl_i))+(sum(sum(sum((b_c-s_c).*dl_c)))));
    if gap<5 && t>100
        break
    end
    gamma=2/(2+t);
    [b_i,b_c,s_i,s_c]=linprog_update( A,cliques,domain_sizes,sample,b_i,b_c,dl_i,dl_c,input,Aeq1,Aeq2,Beq1,Beq2,lb,ub,gamma);

end

predictions=zeros(size(output));
for i=1:length(output)
    [~,predictions(i)]=max(b_i(output(i),:));
end

end

