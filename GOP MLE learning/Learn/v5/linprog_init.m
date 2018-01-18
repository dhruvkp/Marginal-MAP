function [ Aeq1,Aeq2,Beq1,Beq2,lb,ub ] = linprog_init( A,domain_sizes,cliques )
%LINPROG_INIT Summary of this function goes here
%   Detailed explanation goes here
%%%%%  initializing Aeq, beq for linprog - same for all iterations  %%%%%%%

% define order of variables to sum out
number_variables=size(A,1);
number_edges=sum(sum(A))/2;
variables_domain_size=sum(domain_sizes);
edge_domain_size=0;
number_edge_contraints=0;
for c=1:number_edges
    i=cliques(c,1);
    j=cliques(c,2);
    edge_domain_size=edge_domain_size+domain_sizes(i)*domain_sizes(j);
    number_edge_contraints=number_edge_contraints+domain_sizes(i)+domain_sizes(j);
end

all_parameters_size=variables_domain_size+edge_domain_size;

% singleton constraints
Aeq1=zeros(number_variables,all_parameters_size);
Beq1=ones(number_variables,1);
prev_constraints=0;
for i=1:number_variables
    Aeq1(i,(prev_constraints+1):(prev_constraints+domain_sizes(i)))=ones(1,domain_sizes(i));
    prev_constraints=prev_constraints+domain_sizes(i);
end


% sum(cliques)=singleton constraints
Aeq2=zeros(number_edge_contraints,all_parameters_size);
Beq2=zeros(number_edge_contraints,1);
prev_constraints=0;
prev_domain=0;
for c=1:number_edges
    i=cliques(c,1);
    j=cliques(c,2);
    for xi=1:domain_sizes(i)
        %sum xj for xi=xi
        Aeq2(prev_constraints+xi,sum(domain_sizes(1:(i-1)))+xi)=-1;
        for xj=1:domain_sizes(j)
            Aeq2(prev_constraints+xi,(variables_domain_size + prev_domain+domain_sizes(j)*(xi-1)+xj))=1;
        end
    end
    prev_constraints=prev_constraints+domain_sizes(i);
    for xj=1:domain_sizes(j)
        %sum xi for xj=xj
        Aeq2(prev_constraints+xj,sum(domain_sizes(1:(j-1)))+xj)=-1;
        for xi=1:domain_sizes(i)
            Aeq2(prev_constraints+xj,(variables_domain_size + prev_domain+domain_sizes(j)*(xi-1)+xj))=1;
        end
    end
    prev_constraints=prev_constraints+domain_sizes(j);
    prev_domain=prev_domain+domain_sizes(i)*domain_sizes(j);
end


lb=zeros(1,all_parameters_size);
ub=ones(1,all_parameters_size);

end

