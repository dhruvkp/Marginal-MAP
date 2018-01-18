function [b_i,b_c,s_i,s_c,fval]=linprog_update( A,m,domain_size,training_data,b_i,b_c,dl_i,dl_c,fixed,Aeq1,Aeq2,Beq1,Beq2,lb,ub,gamma)
%LINPROG_UPDATE Summary of this function goes here
%   Detailed explanation goes here
% define order of variables to sum out
number_variables=size(A,1);

number_edges=sum(sum(A))/2;


all_parameters_size=number_variables*domain_size+number_edges*(domain_size^2);

    %%% initializing linprog inputs %%%
    s_i=b_i;
    s_c=b_c;
    % extra constraints for setting beliefs of input,output variables to 1
    Aeq3=zeros(length(fixed),all_parameters_size);
    Beq3=ones(length(fixed),1);
    cnt=1;
    for i=fixed
        Aeq3(cnt,domain_size*(i-1)+training_data(m,i))=1;
        cnt=cnt+1;
    end
    
    Aeq=vertcat(Aeq1,Aeq2,Aeq3);
    Beq=vertcat(Beq1,Beq2,Beq3);
    
    % initializing dl -> f in linprog
    dl=zeros(1,number_variables*domain_size+number_edges*(domain_size^2));
    cnt=1;
    for i=1:number_variables
        for xi=1:domain_size
            dl(1,cnt)=dl_i(m,i,xi);
            cnt=cnt+1;
        end
    end
    for c=1:number_edges
        for xi=1:domain_size
            for xj=1:domain_size
                dl(1,cnt)=dl_c(m,c,xi,xj);
                cnt=cnt+1;
            end
        end
    end
    
    % solve linear program for b
    [s_t,fval]=linprog(dl,[],[],Aeq,Beq,lb,ub,[],optimset('Display','none'));
    
    
    cnt=1;
    for i=1:number_variables
        for xi=1:domain_size
            b_i(m,i,xi)=(1-gamma)*b_i(m,i,xi)+gamma*s_t(cnt);
            s_i(m,i,xi)=s_t(cnt);
            cnt=cnt+1;
        end
    end
    for c=1:number_edges
        for xi=1:domain_size
            for xj=1:domain_size
                b_c(m,c,xi,xj)=(1-gamma)*b_c(m,c,xi,xj)+gamma*s_t(cnt);
                s_c(m,c,xi,xj)=s_t(cnt);
                cnt=cnt+1;
            end
        end
    end
end

