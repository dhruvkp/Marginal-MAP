function [b_i,b_c,s_i,s_c,fval]=linprog_update_v2( A,cliques,m,domain_size,training_data,b_i,b_c,dl_i,dl_c,tolearn,number_tolearn_edges,Aeq1,Aeq2,Beq1,Beq2,lb,ub,gamma)
%LINPROG_UPDATE Summary of this function goes here
%   Detailed explanation goes here
% define order of variables to sum out
number_variables=size(A,1);

number_edges=sum(sum(A))/2;


%%% initializing linprog inputs %%%
s_i=b_i;
s_c=b_c;

Aeq=vertcat(Aeq1,Aeq2);
Beq=vertcat(Beq1,Beq2);

% initializing dl -> f in linprog
dl=zeros(1,length(tolearn)*domain_size+number_tolearn_edges*(domain_size^2));
cnt=1;
for i=1:length(tolearn)
    for xi=1:domain_size
        dl(1,cnt)=dl_i(m,tolearn(i),xi);
        cnt=cnt+1;
    end
end
for c=1:number_edges
    i=cliques(c,1);
    j=cliques(c,2);
    [ti,ii]=ismember(i,tolearn);
    [tj,ij]=ismember(j,tolearn);
    if ti && tj
        for xi=1:domain_size
            for xj=1:domain_size
                dl(1,cnt)=dl_c(m,c,xi,xj);
                cnt=cnt+1;
            end
        end
    elseif ti
        for xi=1:domain_size
            dl(1,(ii-1)*domain_size+xi)=dl(1,(ii-1)*domain_size+xi)+dl_c(m,c,xi,training_data(m,j));
        end
    elseif tj
        for xj=1:domain_size
            dl(1,(ij-1)*domain_size+xj)=dl(1,(ij-1)*domain_size+xj)+dl_c(m,c,training_data(m,i),xj);
        end
    end
end

% solve linear program for b
[s_t,fval]=linprog(dl,[],[],Aeq,Beq,lb,ub,[],optimset('Display','none'));


cnt=1;
for i=1:number_variables
    [ti,ii]=ismember(i,tolearn);
    if ti
        for xi=1:domain_size
            b_i(m,i,xi)=(1-gamma)*b_i(m,i,xi)+gamma*s_t(cnt);
            s_i(m,i,xi)=s_t(cnt);
            cnt=cnt+1;
        end
    end
end
for c=1:number_edges
    i=cliques(c,1);
    j=cliques(c,2);
    [ti,ii]=ismember(i,tolearn);
    [tj,ij]=ismember(j,tolearn);
    if ti && tj
        for xi=1:domain_size
            for xj=1:domain_size
                b_c(m,c,xi,xj)=(1-gamma)*b_c(m,c,xi,xj)+gamma*s_t(cnt);
                s_c(m,c,xi,xj)=s_t(cnt);
                cnt=cnt+1;
            end
        end
    elseif ti
        for xi=1:domain_size
            b_c(m,c,xi,training_data(m,j))=(1-gamma)*b_c(m,c,xi,training_data(m,j))+gamma*s_i(m,i,xi);
            s_c(m,c,xi,training_data(m,j))=s_i(m,i,xi);
        end
    elseif tj
        for xj=1:domain_size
            b_c(m,c,training_data(m,i),xj)=(1-gamma)*b_c(m,c,training_data(m,i),xj)+gamma*s_i(m,j,xj);
            s_c(m,c,training_data(m,i),xj)=s_i(m,j,xj);
        end
    end
end
end

