function [ b_i,b_c,bp_i,bp_c ] = init_b(  A,input,output,hidden,domain_sizes,training_data )
%INIT_B initialize beliefs with uniform distribution given training data(evidence)
%   Detailed explanation goes here
M=size(training_data,1);

% define order of variables to sum out
number_variables=size(A,1);

domain_size=max(domain_sizes);
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

b_i=cell(M,1);
for m=1:M
    b_i{m}=sparse(number_variables,domain_size);
end
b_c=cell(M,number_edges);
for m=1:M
    for c=1:number_edges
        b_c{m,c}=sparse(domain_size,domain_size);
    end
end

bp_i=b_i;
bp_c=b_c;

% initialize b_i, bp_i with feasible assignments

for m=1:M
    for i=input
        b_i{m}(i,training_data(m,i))=1;
        bp_i{m}(i,training_data(m,i))=1;
    end
    for i=hidden
        b_i{m}(i,1:domain_sizes(i))=1/domain_sizes(i);
        bp_i{m}(i,1:domain_sizes(i))=1/domain_sizes(i);
    end
    for i=output
        b_i{m}(i,training_data(m,i))=1;
        bp_i{m}(i,1:domain_sizes(i))=1/domain_sizes(i);
    end
end

% initialize b_c, bp_c with feasible assignments

for c=1:number_edges
    i=cliques(c,1);
    j=cliques(c,2);
    if ismember(i,input) && ismember(j,input)
        for m=1:M
            b_c{m,c}(training_data(m,i),training_data(m,j))=1;
            bp_c{m,c}(training_data(m,i),training_data(m,j))=1;
        end
    elseif ismember(i,input) && ~ismember(j,input)
        if ismember(j,output)
            for m=1:M
                b_c{m,c}(training_data(m,i),training_data(m,j))=1;
                bp_c{m,c}(training_data(m,i),1:domain_sizes(j))=1/domain_sizes(j);
            end
        else
            for m=1:M
                b_c{m,c}(training_data(m,i),1:domain_sizes(j))=1/domain_sizes(j);
                bp_c{m,c}(training_data(m,i),1:domain_sizes(j))=1/domain_sizes(j);
            end
        end
    elseif ~ismember(i,input) && ismember(j,input)
        if ismember(i,output)
            for m=1:M
                b_c{m,c}(training_data(m,i),training_data(m,j))=1;
                bp_c{m,c}(1:domain_sizes(i),training_data(m,j))=1/domain_sizes(i);
            end
        else
            for m=1:M
                b_c{m,c}(1:domain_sizes(i),training_data(m,j))=1/domain_sizes(i);
                bp_c{m,c}(1:domain_sizes(i),training_data(m,j))=1/domain_sizes(i);
            end
        end
    else
        for m=1:M
            bp_c{m,c}(1:domain_sizes(i),1:domain_sizes(j))=1/(domain_sizes(i)*domain_sizes(j));
            if ismember(i,output) && ismember(j,output)
                b_c{m,c}(training_data(m,i),training_data(m,j))=1;
            elseif ismember(i,output) && ismember(j,hidden)
                b_c{m,c}(training_data(m,i),1:domain_sizes(j))=1/domain_sizes(j);
            elseif ismember(i,hidden) && ismember(j,output)
                b_c{m,c}(1:domain_sizes(i),training_data(m,j))=1/domain_sizes(i);
            else
                b_c{m,c}(1:domain_sizes(i),1:domain_sizes(j))=1/(domain_sizes(i)*domain_sizes(j));
            end
        end
    end
end

end

