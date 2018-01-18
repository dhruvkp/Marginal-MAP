function [theta_c] = learn( A,input,output,hidden,x_domain,training_data,default_lambda,solver)
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

tolearn=horzcat(hidden,output);

cnt2=0;
for c=1:number_edges
    i=cliques(c,1);
    j=cliques(c,2);
    [ti]=ismember(i,hidden);
    [tj]=ismember(j,hidden);
    if ti && tj
        cnt2=cnt2+1;
    end
end
number_hidden_edges=cnt2;
cnt2=0;
for c=1:number_edges
    i=cliques(c,1);
    j=cliques(c,2);
    [ti]=ismember(i,tolearn);
    [tj]=ismember(j,tolearn);
    if ti && tj
        cnt2=cnt2+1;
    end
end
number_tolearn_edges=cnt2;

%%%% initialize variables to be learned %%%%
 
[ Aeq_1,Aeq_2,Beq_1,Beq_2,lb_1,ub_1 ] = linprog_init( A,domain_size,cliques );
[ Aeq1,Aeq2,Beq1,Beq2,lb,ub ] = linprog_init_v2( A,tolearn,number_tolearn_edges,domain_size,cliques ); 
[ Aeqh1,Aeqh2,Beqh1,Beqh2,lbh,ubh ] = linprog_init_v2( A,hidden,number_hidden_edges,domain_size,cliques ); 

% initializing b, bp
[ b_i,b_c,b_pa_i,bp_i,bp_c,bp_pa_i ] = init_b(  A,input,output,hidden,x_domain,training_data );

% initializing w
[ w_i,w_i_c ] = init_w( A,hidden );

lambda=calc_lambda(A,input,output,hidden,b_i,b_c,w_i,w_i_c,default_lambda,M);

[theta_i,theta_c]=init_theta(A,domain_size);


%%%% learning with frank-wolfe %%%%

s_i=b_i;
s_c=b_c;

sp_i=bp_i;
sp_c=bp_c;

bk_i=b_i;
bk_c=b_c;
bpk_i=bp_i;
bpk_c=bp_c;

bs_i=b_i;
bs_c=b_c;
bps_i=bp_i;
bps_c=bp_c;

t=0;
gap=0; % duality gap
maxim=0;
maxdif=0;

%%% sum of b_i b_c %%%
sum_b_i=zeros(number_variables,domain_size);
sum_b_c=zeros(number_edges,domain_size,domain_size);
for i=1:number_variables
    for xi=1:domain_size
        sum_b_i(i,xi)=sum(b_i(:,i,xi)-bp_i(:,i,xi),1);
    end
end
for c=1:number_edges
    for xi=1:domain_size
        for xj=1:domain_size
            sum_b_c(c,xi,xj)=sum(b_c(:,c,xi,xj)-bp_c(:,c,xi,xj),1);
        end
    end
end

[dl_i,dl_c,dlp_i,dlp_c] = calc_dl_v2( b_i,b_c,bp_i,bp_c,w_i,w_i_c,sum_b_i,sum_b_c,lambda,cliques,order_variables,hidden,output,-1);


while true
    t=t+1;
    %%%% compute duality gap -> check for breaking condition %%%%
        gap=abs(sum(sum(sum((b_i-s_i).*-dl_i)))+sum(sum(sum(sum((b_c-s_c).*-dl_c)))))+abs(sum(sum(sum((bp_i-sp_i).*dlp_i)))+sum(sum(sum(sum((bp_c-sp_c).*dlp_c)))));
    %
            gap
    if t>5000
%         gap
        break;
    end
    %     m=mod(t,M)+1;
    m=randi(M,1);
    gamma=2*M/(2*M+t);
    
    %%% update sum_b_i sum_b_c %%%
    for i=1:number_variables
        for xi=1:domain_size
            sum_b_i(i,xi)=sum_b_i(i,xi)-b_i(m,i,xi)+bp_i(m,i,xi);
        end
    end
    for c=1:number_edges
        for xi=1:domain_size
            for xj=1:domain_size
                sum_b_c(c,xi,xj)=sum_b_c(c,xi,xj)-b_c(m,c,xi,xj)+bp_c(m,c,xi,xj);
            end
        end
    end
    
    
    if strcmp(solver,'linprog2')
        if length(hidden)>0
            [b_i,b_c,s_i,s_c]=linprog_update_v2( A,cliques,m,domain_size,training_data,b_i,b_c,-dl_i,-dl_c,hidden,number_hidden_edges,Aeqh1,Aeqh2,Beqh1,Beqh2,lbh,ubh,gamma);
        end
        [bp_i,bp_c,sp_i,sp_c]=linprog_update_v2( A,cliques,m,domain_size,training_data,bp_i,bp_c,dlp_i,dlp_c,tolearn,number_tolearn_edges,Aeq1,Aeq2,Beq1,Beq2,lb,ub,gamma);
    elseif strcmp(solver,'linprog1')
        if length(hidden)>0
            [b_i,b_c,s_i,s_c]=linprog_update( A,m,domain_size,training_data,b_i,b_c,-dl_i,-dl_c,horzcat(input,output),Aeq_1,Aeq_2,Beq_1,Beq_2,lb_1,ub_1,gamma);
        end
        [bp_i,bp_c,sp_i,sp_c]=linprog_update( A,m,domain_size,training_data,bp_i,bp_c,dlp_i,dlp_c,input,Aeq_1,Aeq_2,Beq_1,Beq_2,lb_1,ub_1,gamma);
    elseif strcmp(solver,'qpbo')
        if length(hidden)>0
            [b_i,b_c,s_i,s_c]=qpbo_update_v2(m,gamma,-dl_i,-dl_c,b_i,b_c,cliques,domain_size,number_variables,number_edges,hidden,number_hidden_edges,training_data(m,:));
        end
        [bp_i,bp_c,sp_i,sp_c]=qpbo_update_v2(m,gamma,dlp_i,dlp_c,bp_i,bp_c,cliques,domain_size,number_variables,number_edges,tolearn,number_tolearn_edges,training_data(m,:));
    elseif strcmp(solver,'compare')
        [dlk_i,dlk_c,dlpk_i,dlpk_c] = calc_dl_v2( bk_i,bk_c,bpk_i,bpk_c,w_i,w_i_c,lambda,cliques,order_variables,hidden,output,m);
        [dls_i,dls_c,dlps_i,dlps_c] = calc_dl_v2( bs_i,bs_c,bps_i,bps_c,w_i,w_i_c,lambda,cliques,order_variables,hidden,output,m );
        [bk_i,bk_c,sk_i,sk_c,fval_linprog]=linprog_update( A,m,domain_size,training_data,bk_i,bk_c,-dlk_i,-dlk_c,horzcat(input,output),Aeq_1,Aeq_2,Beq_1,Beq_2,lb_1,ub_1,gamma);
        [bpk_i,bpk_c,spk_i,spk_c,fvalp_linprog]=linprog_update( A,m,domain_size,training_data,bpk_i,bpk_c,dlpk_i,dlpk_c,input,Aeq_1,Aeq_2,Beq_1,Beq_2,lb_1,ub_1,gamma);
        [bs_i,bs_c,ss_i,ss_c,fvals_linprog]=linprog_update_v2( A,cliques,m,domain_size,training_data,bs_i,bs_c,-dls_i,-dls_c,hidden,number_hidden_edges,Aeqh1,Aeqh2,Beqh1,Beqh2,lbh,ubh,gamma);
        [bps_i,bps_c,sps_i,sps_c,fvalps_linprog]=linprog_update_v2( A,cliques,m,domain_size,training_data,bps_i,bps_c,dlps_i,dlps_c,tolearn,number_tolearn_edges,Aeq1,Aeq2,Beq1,Beq2,lb,ub,gamma);
%         [bs_i,bs_c,ss_i,ss_c,fval_qpbo]=qpbo_update_v2(m,gamma,-dls_i,-dls_c,bs_i,bs_c,cliques,domain_size,number_variables,number_edges,hidden,number_hidden_edges,training_data(m,:));
%         [bps_i,bps_c,sps_i,sps_c,fvalp_qpbo]=qpbo_update_v2(m,gamma,dlps_i,dlps_c,bps_i,bps_c,cliques,domain_size,number_variables,number_edges,tolearn,number_tolearn_edges,training_data(m,:));

        tmp=max(max(abs(sk_i(m,:,:)-ss_i(m,:,:))))+max(max(max(abs(sk_c(m,:,:,:)-ss_c(m,:,:,:)))));
        
        if tmp>0.1
            t1=sum(sum(sk_i(m,:,:).*-dlk_i(m,:,:)))+sum(sum(sum(sk_c(m,:,:,:).*-dlk_c(m,:,:,:))))
            t2=sum(sum(ss_i(m,:,:).*-dls_i(m,:,:)))+sum(sum(sum(ss_c(m,:,:,:).*-dls_c(m,:,:,:))))
            if abs(t1-t2)>1
            end
        end
        
        tmp=max(max(abs(spk_i(m,:,:)-sps_i(m,:,:))))+max(max(max(abs(spk_c(m,:,:,:)-sps_c(m,:,:,:)))));
        
        if tmp>0.1
            t1=sum(sum(spk_i(m,:,:).*dlpk_i(m,:,:)))+sum(sum(sum(spk_c(m,:,:,:).*dlpk_c(m,:,:,:))))
            t2=sum(sum(sps_i(m,:,:).*dlps_i(m,:,:)))+sum(sum(sum(sps_c(m,:,:,:).*dlps_c(m,:,:,:))))
            if abs(t1-t2)>1
            end
        end
        b_i=bk_i;
        b_c=bk_c;
        bp_i=bpk_i;
        bp_c=bpk_c;
        
    end
    
    %%% update sum_b_i sum_b_c %%%
    for i=1:number_variables
        for xi=1:domain_size
            sum_b_i(i,xi)=sum_b_i(i,xi)+b_i(m,i,xi)-bp_i(m,i,xi);
        end
    end
    for c=1:number_edges
        for xi=1:domain_size
            for xj=1:domain_size
                sum_b_c(c,xi,xj)=sum_b_c(c,xi,xj)+b_c(m,c,xi,xj)-bp_c(m,c,xi,xj);
            end
        end
    end
    
    [dln_i,dln_c,dlpn_i,dlpn_c] = calc_dl_v2( b_i,b_c,bp_i,bp_c,w_i,w_i_c,sum_b_i,sum_b_c,lambda,cliques,order_variables,hidden,output,m);
    dl_i(m,:,:)=dln_i(m,:,:);
    dlp_i(m,:,:)=dlpn_i(m,:,:);
    dl_c(m,:,:,:)=dln_c(m,:,:,:);
    dlp_c(m,:,:,:)=dlpn_c(m,:,:,:);
    
end

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