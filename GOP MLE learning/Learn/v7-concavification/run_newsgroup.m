load('C:\Users\drp150030\Downloads\20news_w100.mat');
training_data=full(double(documents'));
% max_line_size=size(training_data,2);
hidden_domain_size=5;
no_nodes=size(A,1);

input=[];
output=1:100;
hidden=101:no_nodes;
domain_sizes=zeros(no_nodes,1);
domain_sizes(output)=2;
domain_sizes(hidden)=2;

number_variables=size(A,1);
order_variables=horzcat(hidden,output,input);
number_edges=sum(sum(A))/2;

sample_size=size(training_data,1);
samples=zeros(sample_size,no_nodes);
samples(:,output)=training_data(1:sample_size,1:100)+1;
samples(:,hidden)=-1;

% for lambda=[0.001,0.01,0.1]
    [theta_i,theta_c_learned]=learn(A,input,output,hidden,domain_sizes,samples,1,'linprog1');
    save('trained_newsgroup.mat');
% end


% 
% roots = t_hat.t0;                  % the list of roots
% children = t_hat.t;                % the list of children for each node
% n = length(children);          % number of nodes
%                                
% parents = zeros(1,n);          % a list of parents
% stack = roots;                 % which nodes need to be processed
% while ~isempty(stack)         
%   node = stack(1);            
%   stack = stack(2:end);        % pop off an element
%   kids = children{node};
%   for kid = kids
%     parents(kid) = node;
%   end
%   stack = [stack, kids];       % push the children of node
% end
% 
% stack=1:100;
% visited=zeros(number_variables,1);
% order_summation=[];
% while ~isempty(stack)
%     node=stack(1);
%     stack=stack(2:end);
%     parent=parents(node);
%     order_summation=[order_summation,node] ;
% end
% 
% logl=0;
% for m=1:M
%     sample=squeeze(samples(m,:));
%     singleton=zeros(number_variables,domain_sizes(1));
%     messages=ones(size(A,1),size(A,1),max(domain_sizes));
% 
%     for i=hidden
%         for xi=1:domain_sizes(i)
%             singleton(i,xi)=exp(theta_c(i,i+1,xi,sample(i+1)));
%         end
%         if i==1
%             singleton(i,xi)=singleton(i,xi)*exp(theta_i(xi));
%         end
%     end
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
%    
% end