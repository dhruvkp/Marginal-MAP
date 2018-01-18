function theta=EM(adjac_mat,theta_size,hidden_indexes,output_indexes,input_indexes,training_data,lambda);

% mark all hidden values as a blank
training_data(:,hidden_indexes)=-1;
% Create a tree out of our graph
input_tree=create_tree(adjac_mat);

number_variables=size(adjac_mat,1);
theta=zeros(theta_size);
x_domain=[1 2];
domain_size=2;

% initialize theta random;y
theta=rand(number_variables,number_variables,domain_size,domain_size)*-1;

for i=1:number_variables
    for j=i:number_variables
        if adjac_mat(i,j)==0
            theta(i,j,:,:)=zeros(domain_size,domain_size);
            theta(j,i,:,:)=theta(i,j,:,:);
        end
        theta(i,j,:,:)=squeeze(theta(j,i,:,:))';
    end
end

    q_of_h=E_step(adjac_mat,input_tree,theta,hidden_indexes,training_data);
    temp_q_of_h=q_of_h+1;
    count=0;
    threshold=ones(size(q_of_h))*.01;
    result1=abs(temp_q_of_h-q_of_h)>threshold;
while  size(find(result1==1),1)>0 && count<100
    result1=abs(temp_q_of_h-q_of_h)>threshold;
    temp_q_of_h=q_of_h;
    count=count+1;
    delta_theta=M_step(adjac_mat,q_of_h,input_tree,theta,hidden_indexes,output_indexes,input_indexes,training_data,lambda);
    theta=theta+0.1*delta_theta;

    q_of_h=E_step(adjac_mat,input_tree,theta,hidden_indexes,training_data);

end
end
