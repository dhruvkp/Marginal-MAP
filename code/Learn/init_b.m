function [ b_i,b_c,b_pa_i,bp_i,bp_c,bp_pa_i ] = init_b(  A,input,output,hidden,x_domain,training_data )
%INIT_B initialize beliefs with uniform distribution given training data(evidence)
%   Detailed explanation goes here
M=size(training_data,1);

% define order of variables to sum out
number_variables=size(A,1);

domain_size=size(x_domain,2);
number_edges=sum(sum(A))/2;

order_variables=horzcat(hidden,output,input);


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


%%%% initialize variables to be learned %%%%


% initialize b with zero

b_i=zeros(M,number_variables,domain_size);
b_c=zeros(M,number_edges,domain_size,domain_size);
b_pa_i=zeros(M,number_edges,domain_size);

bp_i=zeros(M,number_variables,domain_size);
bp_c=zeros(M,number_edges,domain_size,domain_size);
bp_pa_i=zeros(M,number_edges,domain_size);

% initialize b_i, bp_i with feasible assignments

for m=1:M
    for i=input
        b_i(m,i,training_data(m,i))=1;
        bp_i(m,i,training_data(m,i))=1;
    end
    for i=hidden
        b_i(m,i,:)=ones(domain_size,1)/domain_size;
        bp_i(m,i,:)=ones(domain_size,1)/domain_size;
    end
    for i=output
        b_i(m,i,training_data(m,i))=1;
        bp_i(m,i,:)=ones(domain_size,1)/domain_size;
    end
end

% initialize b_c, bp_c with feasible assignments

for c=1:number_edges
    i=cliques(c,1);
    j=cliques(c,2);
    if ismember(i,input) && ismember(j,input)
        for m=1:M
            b_c(m,c,training_data(m,i),training_data(m,j))=1;
            bp_c(m,c,training_data(m,i),training_data(m,j))=1;
        end
    elseif ismember(i,input) && ~ismember(j,input)
        if ismember(j,output)
            for m=1:M
                b_c(m,c,training_data(m,i),training_data(m,j))=1;
                bp_c(m,c,training_data(m,i),:)=ones(domain_size,1)/domain_size;
            end
        else
            for m=1:M
                b_c(m,c,training_data(m,i),:)=ones(domain_size,1)/domain_size;
                bp_c(m,c,training_data(m,i),:)=ones(domain_size,1)/domain_size;
            end
        end
    elseif ~ismember(i,input) && ismember(j,input)
        if ismember(i,output)
            for m=1:M
                b_c(m,c,training_data(m,i),training_data(m,j))=1;
                bp_c(m,c,:,training_data(m,j))=ones(domain_size,1)/domain_size;
            end
        else
            for m=1:M
                b_c(m,c,:,training_data(m,j))=ones(domain_size,1)/domain_size;
                bp_c(m,c,:,training_data(m,j))=ones(domain_size,1)/domain_size;
            end
        end
    else
        for m=1:M
            bp_c(m,c,:,:)=ones(domain_size,domain_size)/(domain_size^2);
            if ismember(i,output) && ismember(j,output)
                b_c(m,c,training_data(m,i),training_data(m,j))=1;
            elseif ismember(i,output) && ismember(j,hidden)
                b_c(m,c,training_data(m,i),:)=ones(domain_size,1)/domain_size;
            elseif ismember(i,hidden) && ismember(j,output)
                b_c(m,c,:,training_data(m,j))=ones(domain_size,1)/domain_size;
            else
                b_c(m,c,:,:)=ones(domain_size,domain_size)/(domain_size^2);
            end
        end
    end
end

% compute b_pa_i, bp_pa_i

for c=1:number_edges
    i=cliques(c,1);
    j=cliques(c,2);
    if find(order_variables==i)<find(order_variables==j) % j is parent of i, j is summed out later than i
        b_pa_i(:,c,:)=sum(b_c(:,c,:,:),3);
        bp_pa_i(:,c,:)=sum(bp_c(:,c,:,:),3);
    else % i is parent of j, i is summed out later than j
        b_pa_i(:,c,:)=sum(b_c(:,c,:,:),4);
        bp_pa_i(:,c,:)=sum(bp_c(:,c,:,:),4);
    end
end
end

