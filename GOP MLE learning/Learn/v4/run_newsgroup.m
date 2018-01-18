% load('C:\Users\drp150030\Downloads\20news_w100.mat');
% training_data=full(double(documents'));
% max_line_size=size(training_data,2);
% hidden_domain_size=5;
% no_nodes=size(A,1);
% 
% input=[];
% output=1:100;
% hidden=101:no_nodes;
% 
% number_variables=size(A,1);
% order_variables=horzcat(hidden,output,input);
% number_edges=sum(sum(A))/2;
% 
% sample_size=size(training_data,1);
% samples=zeros(sample_size,no_nodes);
% samples(:,output)=training_data(1:sample_size,1:100)+1;
% samples(:,hidden)=-1;
% 
%     [theta_i,theta_c_learned]=learn(A,input,output,hidden,[1,2],samples,1,'qpbo');


addpath('../');
roots = t_hat.t0;                  % the list of roots
children = t_hat.t;                % the list of children for each node
n = length(children);          % number of nodes

parents = zeros(1,n);          % a list of parents
stack = roots;                 % which nodes need to be processed
while ~isempty(stack)         
  node = stack(1);            
  stack = stack(2:end);        % pop off an element
  kids = children{node};
  for kid = kids
    parents(kid) = node;
  end
  stack = [stack, kids];       % push the children of node
end
    
max_theta_i=max(max(theta_i))/10;
max_theta_c=max(max(max(max(theta_c_learned))))/10;
theta_i=theta_i-max_theta_i;
theta_c_learned=theta_c_learned-max_theta_c;

[Z,b_i,b_ij]=BP_tree(A,parents,roots,theta_c_learned,theta_i,2*ones(number_variables,1));

theta_i=theta_i+max_theta_i;
theta_c_learned=theta_c_learned+max_theta_c;
Z=Z+max_theta_c+max_theta_i;

logl=0;
for m=1:sample_size
    sample=squeeze(samples(m,:));
    singleton=theta_i;

    for i=output
        j=parents(i);
        for xj=1:2
            singleton(j,xj)=singleton(j,xj)*theta_c_learned(i,j,sample(i),xj);
        end
    end
    A_part=A(hidden,hidden);
    singleton=singleton(hidden,:);
    theta_c=theta_c_learned(hidden,hidden,:,:);
    roots_temp=[];
    for root=roots
        [~,idx]=ismember(root,hidden);
        roots_temp=[roots_temp,idx];
    end
    [~,new_parents]=ismember(parents(hidden),hidden);
    
    
    max_singleton=max(sum(singleton,1));
    max_theta_c=max(max(sum(sum(theta_c,2),1)));
    singleton=singleton-max_singleton;
    theta_c=theta_c-max_theta_c;
    [Z_temp,b_i,b_ij]=BP_tree(A_part,new_parents,roots_temp,theta_c,singleton,2*ones(size(hidden)));
    Z_temp=Z_temp+max_singleton+max_theta_c;
    logl=logl+log(Z_temp/Z);
     log(Z_temp/Z)
end

save('trained_newsgroup_lambda_1.mat');