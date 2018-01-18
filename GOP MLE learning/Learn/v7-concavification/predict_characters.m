function [ predictions ] = predict_characters( A,output,hidden,to_predict,domain_sizes,theta_c,sample )
%PREDICT Summary of this function goes here
%   Detailed explanation goes here
addpath('./permn');
to_predict_combs=permn(1:characters_size,length(to_predict));
% hidden_combs=permn(1:hidden_size,length(hidden));
max_prob=0;
predicitons=zeros(length(to_predict),1);

for i=1:size(to_predict_combs,1)
    sample(to_predict)=to_predict_combs(i,:);
    prob=0;
    for j=1:hidden
        for xj=1:domain_sizes(j)
            prob=prob+exp()
        end
    end
%     for j=1:size(hidden_combs,1)
%         sample(hidden)=hidden_combs(j,:);
%         prob=prob+calc_p_all_x(A,sample,theta_c);
%     end
    if prob>max_prob
        max_prob=prob;
        predicitons(:)=sample(to_predict);
    end
end

% % max-product
% 
% fixed=setdiff(output,to_predict);
% for i=fixed
%     temp=theta_c(i,:,sample(i),:);
%     theta_c(i,:,:,:)=-10^17;
%     theta_c(i,:,sample(i),:)=temp;
%     
%     temp=theta_c(:,i,:,sample(i));
%     theta_c(:,i,:,:)=-10^17;
%     theta_c(:,i,:,sample(i))=temp;
% end
% 
% messages=ones(size(A,1),size(A,1),max(domain_sizes));
% 
% % storing cliques
% number_variables=size(A,1);
% number_edges=sum(sum(A))/2;
% cliques=zeros(number_edges,2);
% cnt=1;
% for i=1:number_variables
%     for j=i+1:number_variables
%         if A(i,j)>0
%             cliques(cnt,:)=[i j];
%             cnt=cnt+1;
%         end
%     end
% end
% 
% % messages- root to leaves
% for c=1:number_edges
%     i=cliques(c,1);
%     j=cliques(c,2);
%     neighbours_i=find(A(i,:));
%     for xj=1:domain_sizes(j)
%         max_m=0;
%         for xi=1:domain_sizes(i)
%             m=theta_c(i,j,xi,xj);
%             for k=neighbours_i
%                 if k==j
%                     continue;
%                 end
%                 m=m*messages(k,i,xi);
%             end
%             if m>max_m
%                 max_m=m;
%             end
%         end
%         messages(i,j,xj)=max_m;
%     end
% end
% 
% % messages- leaves to root
% for c=1:number_edges
%     i=cliques(c,1);
%     j=cliques(c,2);
%     neighbours_j=find(A(j,:));
%     for xi=1:domain_sizes(i)
%         max_m=0;
%         for xj=1:domain_sizes(j)
%             m=theta_c(i,j,xi,xj);
%             for k=neighbours_j
%                 if k==i
%                     continue;
%                 end
%                 m=m*messages(k,j,xj);
%             end
%             if m>max_m
%                 max_m=m;
%             end
%         end
%         messages(j,i,xi)=max_m;
%     end
% end
% 
% predictions=zeros(length(to_predict),1);
% cnt=1;
% for t=to_predict'
%     max_prob=0;
%     max_xt=0;
%     neighbours_t=find(A(t,:));
%     for xt=1:domain_sizes(t)
%         m=1;
%         for n=neighbours_t
%             m=m*messages(n,t,xt);
%         end
%         if m>max_prob
%             max_prob=m;
%             max_xt=xt;
%         end
%     end
%     if max_xt~=20
%     end
%     predictions(cnt)=max_xt;
%     cnt=cnt+1;
% end

end

