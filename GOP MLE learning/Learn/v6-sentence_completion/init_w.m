function [ w_i,w_i_c ] = init_w( A,hidden )
%INIT_W Summary of this function goes here
%   Detailed explanation goes here

number_variables=size(A,1);
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

% initializing constant weights that satisfy constraints

w_i=zeros(number_variables,1);
w_i_c=zeros(number_variables,number_edges);


max_num_cliques_per_variable=0;

for v=hidden
    number_cliques=sum(A(v,hidden));
    if number_cliques>max_num_cliques_per_variable
        max_num_cliques_per_variable=number_cliques;
    end
end
max_num_cliques_per_variable=max_num_cliques_per_variable+1;

for c=1:number_edges
    i=cliques(c,1);
    j=cliques(c,2);
    if ismember(i,hidden) && ismember(j,hidden)
        w_i_c(i,c)=1/max_num_cliques_per_variable;
        w_i_c(j,c)=1/max_num_cliques_per_variable;
    end
end

for v=hidden
    number_cliques=sum(A(v,hidden));
    if max_num_cliques_per_variable~=0
        w_i(v)=1-number_cliques*(1/max_num_cliques_per_variable);
    else
        w_i(v)=1;
    end
end

end

