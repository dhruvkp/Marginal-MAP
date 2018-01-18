function [theta_c] = learn( A,input,output,hidden,x_domain,training_data,default_lambda)
% input: adjacency matrix, indices of input,output and hidden variables,
% domain of input variables, training_data
%
% output: learned parameters for each clique

M=size(training_data,1);

% define order of variables to sum out
number_variables=size(A,1);
order_variables=horzcat(hidden,output,input);

domain_size=size(x_domain,2);
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


%%%% initialize variables to be learned %%%%

% initializing b, bp
[ b_i,b_c,b_pa_i,bp_i,bp_c,bp_pa_i ] = init_b(  A,input,output,hidden,x_domain,training_data );

% initializing w
[ w_i,w_i_c ] = init_w( A,hidden );

lambda=calc_lambda( A,input,output,hidden,b_i,b_c,w_i,w_i_c,default_lambda )

%%%%%  initializing Aeq, beq for linprog - same for all iterations  %%%%%%%


all_parameters_size=number_variables*domain_size+number_edges*(domain_size^2)+number_edges*domain_size;

% singleton constraints
Aeq1=zeros(all_parameters_size,number_variables);
Beq1=ones(1,number_variables);
for i=1:number_variables
    Aeq1((domain_size*(i - 1)+1):(domain_size*i),i)=ones(domain_size,1);
end

% sum(cliques)=singleton constraints
Aeq2=zeros(all_parameters_size,2*domain_size*number_edges);
Beq2=zeros(1,2*domain_size*number_edges);
for c=1:number_edges
    i=cliques(c,1);
    j=cliques(c,2);
    for xi=1:domain_size
        %sum xj for xi=xi
        Aeq2(domain_size*(i-1)+xi,2*domain_size*(c-1)+xi)=-1;
        for xj=1:domain_size
            Aeq2((number_variables*domain_size + (c-1)*(domain_size^2)+domain_size*(xi-1)+xj),2*domain_size*(c-1)+xi)=1;
        end
    end
    for xj=1:domain_size
        %sum xi for xj=xj
        Aeq2(domain_size*(j-1)+xj,2*domain_size*(c-1)+domain_size+xj)=-1;
        for xi=1:domain_size
            Aeq2((number_variables*domain_size + (c-1)*(domain_size^2)+domain_size*(xi-1)+xj),2*domain_size*(c-1)+domain_size+xj)=1;
        end
    end
end

% sum(b_c)=b_pa constraints -> reduces to b_pa=b_i for edge cliques
Aeq3=zeros(all_parameters_size,domain_size*number_edges);
Beq3=zeros(1,domain_size*number_edges);
for c=1:number_edges
    i=cliques(c,1);
    j=cliques(c,2);
    if find(order_variables==i)<find(order_variables==j) % j is parent of i
        for xj=1:domain_size
            Aeq3(domain_size*(j-1)+xj,domain_size*(c-1)+xj)=-1;
            Aeq3(number_variables*domain_size+number_edges*(domain_size^2)+domain_size*(c-1)+xj,domain_size*(c-1)+xj)=1;
        end
    else % i is parent of j
        for xi=1:domain_size
            Aeq3(domain_size*(i-1)+xi,domain_size*(c-1)+xi)=-1;
            Aeq3(number_variables*domain_size+number_edges*(domain_size^2)+domain_size*(c-1)+xi,domain_size*(c-1)+xj)=1;
        end
    end
end

lb=zeros(1,all_parameters_size);
ub=ones(1,all_parameters_size);



%%%% learning with frank-wolfe %%%%


s_t= zeros(M,number_variables*domain_size+number_edges*(domain_size^2)+number_edges*domain_size);
sp_t= zeros(M,number_variables*domain_size+number_edges*(domain_size^2)+number_edges*domain_size);
for m=1:M
    cnt=1;
    for i=1:number_variables
        for xi=1:domain_size
            s_t(m,cnt)=b_i(m,i,xi);
            sp_t(m,cnt)=bp_i(m,i,xi);
            cnt=cnt+1;
        end
    end
    for c=1:number_edges
        for xi=1:domain_size
            for xj=1:domain_size
                s_t(m,cnt)=b_c(m,c,xi,xj);
                sp_t(m,cnt)=bp_c(m,c,xi,xj);
                cnt=cnt+1;
            end
        end
    end
    for c=1:number_edges
        for xi=1:domain_size
            s_t(m,cnt)=b_pa_i(m,c,xi);
            sp_t(m,cnt)=bp_pa_i(m,c,xi);
            cnt=cnt+1;
        end
    end
end

t=0;
gap=0; % duality gap

while true
    t=t+1;
    
    [dl_i,dl_c,dl_pa_i,dlp_i,dlp_c,dlp_pa_i] = calc_dl( b_i,b_c,b_pa_i,bp_i,bp_c,bp_pa_i,w_i,w_i_c,lambda,cliques,order_variables,hidden,output );
    
    
    %%%% compute duality gap -> check for breaking condition %%%%
%     gap=0;
%     for m=1:M
%         cnt=1;
%         for i=1:number_variables
%             for xi=1:domain_size
%                 gap=gap+(b_i(m,i,xi)-s_t(m,cnt))*(-dl_i(m,i,xi))+(bp_i(m,i,xi)-sp_t(m,cnt))*(dlp_i(m,i,xi));
%                 cnt=cnt+1;
%             end
%         end
%         for c=1:number_edges
%             for xi=1:domain_size
%                 for xj=1:domain_size
%                     gap=gap+(b_c(m,c,xi,xj)-s_t(m,cnt))*(-dl_c(m,c,xi,xj))+(bp_c(m,c,xi,xj)-sp_t(m,cnt))*(dlp_c(m,c,xi,xj));
%                     cnt=cnt+1;
%                 end
%             end
%         end
%         for c=1:number_edges
%             for xi=1:domain_size
%                 gap=gap+(b_pa_i(m,c,xi)-s_t(m,cnt))*(-dl_pa_i(m,c,xi))+(bp_pa_i(m,c,xi)-sp_t(m,cnt))*(dlp_pa_i(m,c,xi));
%                 cnt=cnt+1;
%             end
%         end
%     end
    
    if t>3000
        break;
    end
    m=randi(M,1);
    
    %%% initializing linprog inputs %%%
    
    % extra constraints for setting beliefs of input,output variables to 1
    fixed_variables=horzcat(input,output);
    Aeq4=zeros(all_parameters_size,length(fixed_variables));
    Beq4=ones(1,length(fixed_variables));
    cnt=1;
    for i=fixed_variables
        Aeq4(domain_size*(i-1)+training_data(m,i),cnt)=1;
        cnt=cnt+1;
    end
    
    
    Aeq=horzcat(Aeq1,Aeq2,Aeq3,Aeq4);
    Beq=horzcat(Beq1,Beq2,Beq3,Beq4);
    Aeq=Aeq';
    Beq=Beq';
    
    % initializing dl -> f in linprog
    dl=zeros(1,number_variables*domain_size+number_edges*(domain_size^2)+number_edges*domain_size);
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
    for c=1:number_edges
        for xi=1:domain_size
            dl(1,cnt)=dl_pa_i(m,c,xi);
            cnt=cnt+1;
        end
    end
    
    % solve linear program for b
    s_t(m,:)=linprog(-dl,[],[],Aeq,Beq,lb,ub,[],optimset('Display','none'));
    
    gamma=2*M/(2*M+t);
    
    cnt=1;
    for i=1:number_variables
        for xi=1:domain_size
            b_i(m,i,xi)=(1-gamma)*b_i(m,i,xi)+gamma*s_t(m,cnt);
            cnt=cnt+1;
        end
    end
    for c=1:number_edges
        for xi=1:domain_size
            for xj=1:domain_size
                b_c(m,c,xi,xj)=(1-gamma)*b_c(m,c,xi,xj)+gamma*s_t(m,cnt);
                cnt=cnt+1;
            end
        end
    end
    for c=1:number_edges
        for xi=1:domain_size
            b_pa_i(m,c,xi)=(1-gamma)*b_pa_i(m,c,xi)+gamma*s_t(m,cnt);
            cnt=cnt+1;
        end
    end
    
    % extra constraints for setting beliefs of input variables to 1
    Aeq4=zeros(all_parameters_size,length(input));
    Beq4=ones(1,length(input));
    cnt=1;
    for i=input
        Aeq4(domain_size*(i-1)+training_data(m,i),cnt)=1;
        cnt=cnt+1;
    end
    
    
    Aeq=horzcat(Aeq1,Aeq2,Aeq3,Aeq4);
    Beq=horzcat(Beq1,Beq2,Beq3,Beq4);
    Aeq=Aeq';
    Beq=Beq';
    
    % initializing dl -> f in linprog
    dlp=zeros(1,number_variables*domain_size+number_edges*(domain_size^2)+number_edges*domain_size);
    cnt=1;
    for i=1:number_variables
        for xi=1:domain_size
            dlp(1,cnt)=dlp_i(m,i,xi);
            cnt=cnt+1;
        end
    end
    for c=1:number_edges
        for xi=1:domain_size
            for xj=1:domain_size
                dlp(1,cnt)=dlp_c(m,c,xi,xj);
                cnt=cnt+1;
            end
        end
    end
    for c=1:number_edges
        for xi=1:domain_size
            dlp(1,cnt)=dlp_pa_i(m,c,xi);
            cnt=cnt+1;
        end
    end
    
    % solve linear program for bp
    sp_t(m,:)=linprog(dlp,[],[],Aeq,Beq,lb,ub,[],optimset('Display','none'));
    
    gamma=2*M/(2*M+t);
    
    cnt=1;
    for i=1:number_variables
        for xi=1:domain_size
            bp_i(m,i,xi)=(1-gamma)*bp_i(m,i,xi)+gamma*sp_t(m,cnt);
            cnt=cnt+1;
        end
    end
    for c=1:number_edges
        for xi=1:domain_size
            for xj=1:domain_size
                bp_c(m,c,xi,xj)=(1-gamma)*bp_c(m,c,xi,xj)+gamma*sp_t(m,cnt);
                cnt=cnt+1;
            end
        end
    end
    for c=1:number_edges
        for xi=1:domain_size
            bp_pa_i(m,c,xi)=(1-gamma)*bp_pa_i(m,c,xi)+gamma*sp_t(m,cnt);
            cnt=cnt+1;
        end
    end
    
end

M
t
theta_i=zeros(number_variables,domain_size);
theta_c=zeros(number_variables,number_variables,domain_size,domain_size);

for i=1:number_variables
    for xi=1:domain_size
        theta_i(i,xi)=(sum(b_i(:,i,xi))-sum(bp_i(:,i,xi)))/lambda;
    end
end
for c=1:number_edges
    i=cliques(c,1);
    j=cliques(c,2);
    for xi=1:domain_size
        for xj=1:domain_size
            theta_c(i,j,xi,xj)=(sum(b_c(:,c,xi,xj))-sum(bp_c(:,c,xi,xj)))/lambda;
            theta_c(j,i,xj,xi)=theta_c(i,j,xi,xj);
        end
    end
end
end