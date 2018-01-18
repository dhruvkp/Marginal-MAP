function [theta_i,theta_c_n] = learn( A,input,output,hidden,domain_sizes,training_data,default_lambda,solver)
% input: adjacency matrix, indices of input,output and hidden variables,
% domain of input variables, training_data
%
% output: learned parameters for each clique

M=size(training_data,1);

% define order of variables to sum out
number_variables=size(A,1);
order_variables=horzcat(hidden,output,input);

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

transition_cliques=[];
emission_cliques=[];
number_hidden_edges=0;
number_tolearn_edges=number_edges;
for c=1:number_edges
    i=cliques(c,1);
    j=cliques(c,2);
    [ti]=ismember(i,hidden);
    [tj]=ismember(j,hidden);
    if ti && tj
        number_hidden_edges=number_hidden_edges+1;
        transition_cliques=[transition_cliques;c];
    elseif ti && ismember(j,output)
        emission_cliques=[emission_cliques;c];
    end
end


%%%% initialize variables to be learned %%%%

[ Aeq1,Aeq2,Beq1,Beq2,lb,ub ] = linprog_init( A,domain_sizes,cliques );

% initializing b, bp
[b_i,b_c,bp_i,bp_c] = init_b(  A,input,output,hidden,domain_sizes,training_data );

% initializing w
[w_i,w_i_c] = init_w( A,hidden );

lambda=calc_lambda(A,input,output,hidden,b_i,b_c,w_i,w_i_c,default_lambda,M,domain_sizes);
% lambda=default_lambda;

[ theta_i,theta_c ] = init_theta( A,domain_sizes );


%%%% learning with frank-wolfe %%%%

s_i=b_i;
s_c=b_c;

sp_i=bp_i;
sp_c=bp_c;


t=0;
gap=0; % duality gap
maxim=0;
maxdif=0;

max_domain_size=max(domain_sizes);

%%% sum of b_i b_c %%%
sum_b_i=zeros(number_variables,max_domain_size);
sum_b_c=cell(number_edges,1);
for c=1:number_edges
    sum_b_c{c}=sparse(max_domain_size,max_domain_size);
end

for i=1:number_variables
    for xi=1:domain_sizes(i)
        for m=1:M
            sum_b_i(i,xi)=sum_b_i(i,xi)+b_i{m}(i,xi);
        end
    end
end
for c=1:number_edges
    i=cliques(c,1);
    j=cliques(c,2);
    for xi=1:domain_sizes(i)
        for xj=1:domain_sizes(j)
            for m=1:M
                sum_b_c{c}(xi,xj)=sum_b_c{c}(xi,xj)+b_c{m,c}(xi,xj)-bp_c{m,c}(xi,xj);
            end
        end
    end
end

[dl_i,dl_c,dlp_i,dlp_c,dl_theta_i,dl_theta_c] = calc_dl_v2( b_i,b_c,bp_i,bp_c,theta_i,theta_c,w_i,w_i_c,sum_b_i,sum_b_c,lambda,cliques,order_variables,hidden,output,domain_sizes,-1);

step_theta=1/max([max(abs(dl_theta_c(:))),max(abs(dl_theta_i(:)))]);


while true
    t=t+1;
    %%%% compute duality gap -> check for breaking condition %%%%
    %     gap=abs(sum(sum(sum((b_i-s_i).*-dl_i)))+sum(sum(sum(sum((b_c-s_c).*-dl_c)))))+abs(sum(sum(sum((bp_i-sp_i).*dlp_i)))+sum(sum(sum(sum((bp_c-sp_c).*dlp_c)))));
    %
    %     if isnan(gap)
    %     end
    %     theta_gap=sum(sum(abs(dl_theta_i)))+sum(sum(sum(abs(dl_theta_c))));
    %     if t>10000
    %         gap
    %         theta_gap
    %         break;
    %     end
    if t>81000
        %         gap
        %         theta_gap
        break
    end
    %     m=mod(t,M)+1;
    m=randi(M,1);
    gamma=2*M/(2*M+t);
    
    
    %%% update sum_b_i sum_b_c %%%
    for i=1:number_variables
        for xi=1:domain_sizes(i)
            sum_b_i(i,xi)=sum_b_i(i,xi)-b_i{m}(i,xi)+bp_i{m}(i,xi);
        end
    end
    for c=1:number_edges
        i=cliques(c,1);
        j=cliques(c,2);
        for xi=1:domain_sizes(i)
            for xj=1:domain_sizes(j)
                sum_b_c{c}(xi,xj)=sum_b_c{c}(xi,xj)-b_c{m,c}(xi,xj)+bp_c{m,c}(xi,xj);
            end
        end
    end
    
    
    if ~isempty(hidden)
        [b_i,b_c,s_i,s_c]=linprog_update( A,m,cliques,domain_sizes,training_data(m,:),b_i,b_c,dl_i,dl_c,horzcat(input,output),Aeq1,Aeq2,Beq1,Beq2,lb,ub,gamma);
    end
    [bp_i,bp_c,sp_i,sp_c]=linprog_update( A,m,cliques,domain_sizes,training_data(m,:),bp_i,bp_c,dlp_i,dlp_c,input,Aeq1,Aeq2,Beq1,Beq2,lb,ub,gamma);
    
    
    theta_i_n=theta_i+step_theta*dl_theta_i;
    theta_c=theta_c+step_theta*dl_theta_c;
%     
%     %%% projected gradient ascent- average theta that need to be same
%     for xi=1:domain_sizes(1)
%         theta_i_n(hidden,xi)=sum(theta_i_n(hidden,xi),1)/length(hidden);
%     end
%     for xi=1:domain_sizes(2)
%         theta_i_n(output,xi)=sum(theta_i_n(output,xi),1)/length(output);
%     end
% %     theta_i(1,:)=theta_i_n(1,:);
%     theta_i=theta_i_n;
% 
%     for xi=1:domain_sizes(1)
%         for xj=1:domain_sizes(1)
%             theta_c(transition_cliques,xi,xj)=sum(theta_c(transition_cliques,xi,xj),1)/number_hidden_edges;
%         end
%     end
%     for xi=1:domain_sizes(2)
%         for xj=1:domain_sizes(2)
%             theta_c(emission_cliques,xi,xj)=sum(theta_c(emission_cliques,xi,xj),1)/(number_edges-number_hidden_edges);
%         end
%     end
%     
    
    
    
    %%% update sum_b_i sum_b_c %%%
    for i=1:number_variables
        for xi=1:domain_sizes(i)
            sum_b_i(i,xi)=sum_b_i(i,xi)+b_i{m}(i,xi)-bp_i{m}(i,xi);
        end
    end
    for c=1:number_edges
        i=cliques(c,1);
        j=cliques(c,2);
        for xi=1:domain_sizes(i)
            for xj=1:domain_sizes(j)
                sum_b_c{c}(xi,xj)=sum_b_c{c}(xi,xj)+b_c{m,c}(xi,xj)-bp_c{m,c}(xi,xj);
            end
        end
    end
    
    [dln_i,dln_c,dlpn_i,dlpn_c,dl_theta_i,dl_theta_c] = calc_dl_v2( b_i,b_c,bp_i,bp_c,theta_i,theta_c,w_i,w_i_c,sum_b_i,sum_b_c,lambda,cliques,order_variables,hidden,output,domain_sizes,m);
    dl_i{m}(:,:)=dln_i{m}(:,:);
    dlp_i{m}(:,:)=dlpn_i{m}(:,:);
    for c=1:number_edges
        dl_c{m,c}(:,:)=dln_c{m,c}(:,:);
        dlp_c{m,c}(:,:)=dlpn_c{m,c}(:,:);
    end
end

% theta_i=zeros(number_variables,domain_size);
theta_c_n=zeros(number_variables,number_variables,max_domain_size,max_domain_size);
theta_i_first=theta_i(1,:);
% for i=1:number_variables
%     for xi=1:domain_size
%         theta_i(i,xi)=(sum(b_i(:,i,xi))-sum(bp_i(:,i,xi)))/lambda;
%     end
% end
for c=1:number_edges
    i=cliques(c,1);
    j=cliques(c,2);
    for xi=1:domain_sizes(i)
        for xj=1:domain_sizes(j)
            theta_c_n(i,j,xi,xj)=theta_c(c,xi,xj);
            theta_c_n(j,i,xj,xi)=theta_c(c,xi,xj);
        end
    end
end
end