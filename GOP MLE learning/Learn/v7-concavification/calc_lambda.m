function [ lambda ] = calc_lambda( A,input,output,hidden,b_i,b_c,w_i,w_i_c,default_lambda,M,domain_sizes )
%CALC_LAMBDA calculates regularizer lambda value

number_variables=size(A,2);
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

% finding lambda
lambda=0;
for i=hidden
    max_b_i=0;
    for m=1:M
        for xi=1:domain_sizes(i)
            if b_i{m}(i,xi)>max_b_i
                max_b_i=b_i{m}(i,xi);
            end
        end
    end
    tmp=(M+1)*max_b_i/(2*w_i(i,1));
    if tmp>lambda && ~isinf(tmp)
        lambda=tmp;
    end
end
for c=1:number_edges
    i=cliques(c,1);
    j=cliques(c,2);
    tmp=0;
    max_b_c=0;
    for m=1:M
        for xi=1:domain_sizes(i)
            for xj=1:domain_sizes(j)
                if b_c{m,c}(xi,xj)>max_b_c
                    max_b_c=b_c{m,c}(xi,xj);
                end
            end
        end
    end
    if find(order_variables==i)<find(order_variables==j) && ismember(i,hidden)
        tmp=(M+1)*max_b_c/(2*w_i_c(i,c));
    elseif ismember(j,hidden)
        tmp=(M+1)*max_b_c/(2*w_i_c(j,c));
    end
    if tmp>lambda && ~isinf(tmp)
        lambda=tmp;
    end
end
if lambda<default_lambda
    lambda=default_lambda;
else
    lambda=lambda+0.5;
end


end

