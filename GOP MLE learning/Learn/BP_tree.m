function [Z,b_i,b_ij]=BP_tree(A,pred,roots,theta_c,theta_i,domain_sizes)

messages=ones(size(A,1),size(A,1),max(domain_sizes));
disc=[];
for root=roots
    [disc1, pred1, closed] = graphtraverse(sparse(A),root,'Directed',false,'Method','BFS');
    disc=[disc,disc1];
end
disc=fliplr(disc);
pred=pred(disc);

number_variables=size(A,1);
number_edges=sum(sum(A))/2;
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

% messages- leaves to root
for cnt=1:length(disc)
    i=disc(cnt);
    j=pred(cnt);
    if isnan(j) || j==0
        continue
    end
    neighbours_i=find(A(i,:));
    for xj=1:domain_sizes(j)
        total_m=0;
        for xi=1:domain_sizes(i)
            m=exp(theta_i(i,xi))*exp(theta_c(i,j,xi,xj));
            for k=neighbours_i
                if k==j
                    continue;
                end
                m=m*messages(k,i,xi);
            end
            total_m=total_m+m;
        end
        messages(i,j,xj)=total_m;
    end
    if sum(messages(i,j,:))==0
        messages(i,j,:)=0;
    else
        messages(i,j,:)=messages(i,j,:)/sum(messages(i,j,:));
    end
end

disc=fliplr(disc);
pred=fliplr(pred);

% messages- root to leaves
for cnt=1:length(disc)
    i=disc(cnt);
    parent_i=pred(cnt);
    neighbours_i=find(A(i,:));
    for j=neighbours_i
        if j==parent_i
            continue
        end
        for xj=1:domain_sizes(j)
            total_m=0;
            for xi=1:domain_sizes(i)
                m=exp(theta_i(i,xi))*exp(theta_c(i,j,xi,xj));
                for k=neighbours_i
                    if k==j
                        continue;
                    end
                    m=m*messages(k,i,xi);
                end
                total_m=total_m+m;
            end
            messages(i,j,xj)=total_m;
        end
        if sum(messages(i,j,:))==0
            messages(i,j,:)=0;
        else
            messages(i,j,:)=messages(i,j,:)/sum(messages(i,j,:));
        end
    end
end

b_i=zeros(number_variables,max(domain_sizes));

for i=1:number_variables
    for xi=1:domain_sizes(i)
        m=exp(theta_i(i,xi));
        neighbours_i=find(A(i,:));
        for k=neighbours_i
            m=m*messages(k,i,xi);
        end
        b_i(i,xi)=m;
    end
    b_i(i,:)=b_i(i,:)/sum(b_i(i,:));
end


b_ij=zeros(number_variables,number_variables,max(domain_sizes),max(domain_sizes));

for c=1:number_edges
    i=cliques(c,1);
    j=cliques(c,2);
    neighbours_i=find(A(i,:));
    neighbours_j=find(A(j,:));
    for xi=1:domain_sizes(i)
        for xj=1:domain_sizes(j)
            m=exp(theta_i(i,xi))*exp(theta_i(j,xj))*exp(theta_c(i,j,xi,xj));
            for k=neighbours_i
                if k==j
                    continue
                end
                m=m*messages(k,i,xi);
            end
            for k=neighbours_j
                if k==i
                    continue
                end
                m=m*messages(k,j,xj);
            end
            b_ij(i,j,xi,xj)=m;
        end
    end
    b_ij(i,j,:,:)=b_ij(i,j,:,:)/sum(sum(b_ij(i,j,:,:)));
end

Z=0;


for i=1:number_variables
    for xi=1:domain_sizes(i)
        if b_i(i,xi)~=0
            Z=Z-b_i(i,xi)*log(b_i(i,xi));
        end
        Z=Z+b_i(i,xi)*theta_i(i,xi);
    end
end
for c=1:number_edges
    i=cliques(c,1);
    j=cliques(c,2);
    for xi=1:domain_sizes(i)
        for xj=1:domain_sizes(j)
            if b_ij(i,j,xi,xj)~=0
                Z=Z-b_ij(i,j,xi,xj)*log(b_ij(i,j,xi,xj)/(b_i(i,xi)*b_i(j,xj)));
            end
            Z=Z+b_ij(i,j,xi,xj)*theta_c(i,j,xi,xj);
        end
    end
end

end