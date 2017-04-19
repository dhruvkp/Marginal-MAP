function [ lambda ] = calc_lambda( A,input,output,hidden,b_i,b_c,w_i,w_i_c )
%CALC_LAMBDA calculates regularizer lambda value

number_variables=size(b_i,2);
number_edges=size(b_c,2);

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
for i=1:number_variables
    tmp=0;
    if w_i(i,1)~=0
        tmp=max(max(b_i(:,i,:)))/w_i(i,1);
    end
    if tmp>lambda
        lambda=tmp;
    end
end
for c=1:number_edges
    i=cliques(c,1);
    j=cliques(c,2);
    tmp=0;
    if find(order_variables==i)<find(order_variables==j) && w_i_c(i,c)~=0
        tmp=max(max(max(b_c(:,c,:,:))))/w_i_c(i,c);
    elseif w_i_c(j,c)~=0
        tmp=max(max(max(b_c(:,c,:,:))))/w_i_c(j,c);
    end
    if tmp>lambda
        lambda=tmp;
    end
end


lambda=lambda+0.1;

end

