function [ predictions ] = bp_predict_characters( A,output,input,hidden,to_predict,domain_sizes,theta_i,theta_c,sample )
%PREDICT Summary of this function goes here
%   Detailed explanation goes here
% define order of variables to sum out
characters_size=domain_sizes(2);
to_predict_combs=permn(1:characters_size,length(to_predict));
fixed=[output,to_predict];

% sum-product -> compute Z
theta_offset=max(theta_c(:))/2;
% theta_c=theta_c-theta_offset;

messages=ones(size(A,1),size(A,1),max(domain_sizes));

% storing cliques
number_variables=size(A,1);
number_edges=sum(sum(A))/2;
cliques=zeros(number_edges,2);
cnt=1;
for i=1:number_variables
    j=i+2;
    if j<number_variables && A(i,j)>0
        cliques(cnt,:)=[i j];
        cnt=cnt+1;
    end
end


%%% sum out hidden to compute probability of observed- variable elimination
cnt=1;
max_logp=0;
predictions=to_predict_combs(1,:);
for to_predict_values=to_predict_combs'
    cnt=cnt+1;
    logp=0;
    sample(to_predict)=to_predict_values;
    singleton=zeros(number_variables,domain_sizes(1));
    for i=hidden
        for xi=1:domain_sizes(i)
            singleton(i,xi)=exp(theta_c(i,i+1,xi,sample(i+1)));
        end
        if i==1
            singleton(i,xi)=singleton(i,xi)*exp(theta_i(xi));
        end
    end
    
%     
%     % messages- leaves to root
%     for c=number_edges:-1:1
%         i=cliques(c,1);
%         j=cliques(c,2);
%         neighbours_j=find(A(j,:));
%         for xi=1:domain_sizes(i)
%             total_m=0;
%             for xj=1:domain_sizes(j)
%                 m=exp(theta_c(i,j,xi,xj));
%                 for k=neighbours_j
%                     if k==i
%                         continue;
%                     end
%                     m=m*messages(k,j,xj);
%                 end
%                 total_m=total_m+m;
%             end
%             messages(j,i,xi)=total_m;
%         end
%     end
%     
%     
%     
%     % messages- root to leaves
%     for c=1:number_edges
%         i=cliques(c,1);
%         j=cliques(c,2);
%         neighbours_i=find(A(i,:));
%         for xj=1:domain_sizes(j)
%             total_m=0;
%             for xi=1:domain_sizes(i)
%                 m=exp(theta_c(i,j,xi,xj));
%                 for k=neighbours_i
%                     if k==j
%                         continue;
%                     end
%                     m=m*messages(k,i,xi);
%                 end
%                 total_m=total_m+m;
%             end
%             messages(i,j,xj)=total_m;
%         end
%     end
%     
%     marginals=zeros(number_variables,max(domain_sizes));
%     
%     for i=1:number_variables
%         for xi=1:domain_sizes(i)
%             m=1;
%             neighbours_i=find(A(i,:));
%             for k=neighbours_i
%                 m=m*messages(k,i,xi);
%             end
%             marginals(i,xi)=m;
%         end
%     end
%     
%     logz=log(sum(marginals(1,:)))+theta_offset;


    
    %%% variable elimination
    for i=hidden
        if i==number_variables-1
            for xi=1:domain_sizes(i)
                logp=logp+singleton(i,xi);
            end
            break
        end
        j=i+2;
        for xj=1:domain_sizes(j)
            total=0;
            for xi=1:domain_sizes(i)
                total=total+singleton(i,xi)*exp(theta_c(i,j,xi,xj));
            end
            singleton(j,xj)=singleton(j,xj)*total;
        end
    end
    %     logp=logp-logz;

    
    if logp>max_logp
        max_logp=logp;
        predictions=to_predict_values;
    end
end

% addpath('../../../EM/V5');
% sample(to_predict)=randi(domain_sizes(2),length(to_predict),1);
% sample(hidden)=-1;
% fixed=[output,to_predict];
% [ input_tree , ~ ]=create_tree1(A ,domain_sizes);
% p_of_observed=inference1(input_tree,theta_c,sample);


end