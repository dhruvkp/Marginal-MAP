function [ Aeq1,Aeq2,Beq1,Beq2,lb,ub ] = linprog_init( A,domain_size,cliques )
%LINPROG_INIT Summary of this function goes here
%   Detailed explanation goes here
%%%%%  initializing Aeq, beq for linprog - same for all iterations  %%%%%%%

% define order of variables to sum out
number_variables=size(A,1);
number_edges=sum(sum(A))/2;


all_parameters_size=number_variables*domain_size+number_edges*(domain_size^2);

% singleton constraints
Aeq1=zeros(number_variables,all_parameters_size);
Beq1=ones(number_variables,1);
for i=1:number_variables
    Aeq1(i,(domain_size*(i - 1)+1):(domain_size*i))=ones(1,domain_size);
end

% sum(cliques)=singleton constraints
Aeq2=zeros(2*domain_size*number_edges,all_parameters_size);
Beq2=zeros(2*domain_size*number_edges,1);
for c=1:number_edges
    i=cliques(c,1);
    j=cliques(c,2);
    for xi=1:domain_size
        %sum xj for xi=xi
        Aeq2(2*domain_size*(c-1)+xi,domain_size*(i-1)+xi)=-1;
        for xj=1:domain_size
            Aeq2(2*domain_size*(c-1)+xi,(number_variables*domain_size + (c-1)*(domain_size^2)+domain_size*(xi-1)+xj))=1;
        end
    end
    for xj=1:domain_size
        %sum xi for xj=xj
        Aeq2(2*domain_size*(c-1)+domain_size+xj,domain_size*(j-1)+xj)=-1;
        for xi=1:domain_size
            Aeq2(2*domain_size*(c-1)+domain_size+xj,(number_variables*domain_size + (c-1)*(domain_size^2)+domain_size*(xi-1)+xj))=1;
        end
    end
end
lb=zeros(1,all_parameters_size);
ub=ones(1,all_parameters_size);

end

