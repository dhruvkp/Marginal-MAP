function [b_i,b_c,s_i,s_c,fval]=linprog_update( A,m,cliques,domain_sizes,training_data,b_i,b_c,dl_i,dl_c,fixed,Aeq1,Aeq2,Beq1,Beq2,lb,ub,gamma)
%LINPROG_UPDATE Summary of this function goes here
%   Detailed explanation goes here
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

    %%% initializing linprog inputs %%%
    s_i=b_i;
    s_c=b_c;
    % extra constraints for setting beliefs of input,output variables to 1
    Aeq3=zeros(length(fixed),all_parameters_size);
    Beq3=ones(length(fixed),1);
    cnt=1;
    for i=fixed
        Aeq3(cnt,sum(domain_sizes(1:i-1))+training_data(m,i))=1;
        cnt=cnt+1;
    end
    
    Aeq=vertcat(Aeq1,Aeq2,Aeq3);
    Beq=vertcat(Beq1,Beq2,Beq3);
    
    % initializing dl -> f in linprog
    dl=zeros(1,all_parameters_size);
    cnt=1;
    for i=1:number_variables
        for xi=1:domain_sizes(i)
            dl(1,cnt)=dl_i(m,i,xi);
            cnt=cnt+1;
        end
    end
    for c=1:number_edges
        i=cliques(c,1);
        j=cliques(c,2);
        for xi=1:domain_sizes(i)
            for xj=1:domain_sizes(j)
                dl(1,cnt)=dl_c(m,c,xi,xj);
                cnt=cnt+1;
            end
        end
    end
    
    % solve linear program for b
    [s_t,fval]=linprog(dl,[],[],Aeq,Beq,lb,ub,[],optimset('Display','none'));
    
    
    cnt=1;
    for i=1:number_variables
        for xi=1:domain_sizes(i)
            b_i(m,i,xi)=(1-gamma)*b_i(m,i,xi)+gamma*s_t(cnt);
            s_i(m,i,xi)=s_t(cnt);
            cnt=cnt+1;
        end
    end
    for c=1:number_edges
        i=cliques(c,1);
        j=cliques(c,2);
        for xi=1:domain_sizes(i)
            for xj=1:domain_sizes(j)
                b_c(m,c,xi,xj)=(1-gamma)*b_c(m,c,xi,xj)+gamma*s_t(cnt);
                s_c(m,c,xi,xj)=s_t(cnt);
                cnt=cnt+1;
            end
        end
    end
end

