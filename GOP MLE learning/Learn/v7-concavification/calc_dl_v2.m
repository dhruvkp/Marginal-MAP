function [ dl_i,dl_c,dlp_i,dlp_c,dl_theta_i,dl_theta_c ] = calc_dl_v2( b_i,b_c,bp_i,bp_c,theta_i,theta_c,degree,sum_b_i,sum_b_c,lambda,cliques,order_variables,hidden,output,domain_sizes,m )
%CALC_DL calculates delta-L
%   Detailed explanation goes here

dl_i=b_i;
dl_c=b_c;

dlp_i=b_i;
dlp_c=b_c;

dl_theta_i=zeros(size(theta_i));
dl_theta_c=zeros(size(theta_c));


M=size(b_i,1);
number_variables=size(b_i{1},1);
number_edges=size(b_c,2);

if m==-1
    for m=1:M
        for i=1:number_variables
            for xi=1:domain_sizes(i)
                dl_i{m}(i,xi)=theta_i(i,xi)+(1-2*b_i{m}(i,xi))/(2*lambda/M);
                if b_i{m}(i,xi)~=0
                    dl_i{m}(i,xi)=dl_i{m}(i,xi) - (1- degree(i))*(1+log(b_i{m}(i,xi)));
                end
                dlp_i{m}(i,xi)=-theta_i(i,xi);
                if bp_i{m}(i,xi)~=0
                    dlp_i{m}(i,xi)=dlp_i{m}(i,xi)+(1- degree(i))*(1+log(bp_i{m}(i,xi)));
                end
            end
        end
        
        for c=1:number_edges
            i=cliques(c,1);
            j=cliques(c,2);
            for xi=1:domain_sizes(i)
                for xj=1:domain_sizes(j)
                    dl_c{m,c}(xi,xj)=theta_c(c,xi,xj)+(1-2*b_c{m,c}(xi,xj))/(2*lambda/M);
                    dlp_c{m,c}(xi,xj)=-theta_c(c,xi,xj);
                    if b_c{m,c}(xi,xj)~=0
                        dl_c{m,c}(xi,xj)=dl_c{m,c}(xi,xj) - 1-  log(b_c{m,c}(xi,xj));
                    end
                    if bp_c{m,c}(xi,xj)~=0
                        dlp_c{m,c}(xi,xj)=dlp_c{m,c}(xi,xj)+ 1+ log(bp_c{m,c}(xi,xj));
                    end
                    
                end
            end
        end
        for c=1:number_edges
            dl_c{m,c}(:,:)=-dl_c{m,c}(:,:);
        end
        dl_i{m}(:,:)=-dl_i{m}(:,:);
    end
else
    for i=1:number_variables
        for xi=1:domain_sizes(i)
            dl_i{m}(i,xi)=theta_i(i,xi)+(1-2*b_i{m}(i,xi))/(2*lambda/M);
            if b_i{m}(i,xi)~=0
                dl_i{m}(i,xi)=dl_i{m}(i,xi) - (1- degree(i))*(1+log(b_i{m}(i,xi)));
            end
            dlp_i{m}(i,xi)=-theta_i(i,xi);
            if bp_i{m}(i,xi)~=0
                dlp_i{m}(i,xi)=dlp_i{m}(i,xi)+(1- degree(i))*(1+log(bp_i{m}(i,xi)));
            end
        end
    end
    
    for c=1:number_edges
        i=cliques(c,1);
        j=cliques(c,2);
        for xi=1:domain_sizes(i)
            for xj=1:domain_sizes(j)
                dl_c{m,c}(xi,xj)=theta_c(c,xi,xj)+(1-2*b_c{m,c}(xi,xj))/(2*lambda/M);
                dlp_c{m,c}(xi,xj)=-theta_c(c,xi,xj);
                if b_c{m,c}(xi,xj)~=0
                    dl_c{m,c}(xi,xj)=dl_c{m,c}(xi,xj) - 1-  log(b_c{m,c}(xi,xj));
                end
                if bp_c{m,c}(xi,xj)~=0
                    dlp_c{m,c}(xi,xj)=dlp_c{m,c}(xi,xj)+ 1+ log(bp_c{m,c}(xi,xj));
                end
                
            end
        end
    end
    for c=1:number_edges
        dl_c{m,c}(:,:)=-dl_c{m,c}(:,:);
    end
    dl_i{m}(:,:)=-dl_i{m}(:,:);
end

for i=1:number_variables
    for xi=1:domain_sizes(i)
        dl_theta_i(i,xi)=sum_b_i(i,xi)-lambda*theta_i(i,xi);
    end
end
for c=1:number_edges
    i=cliques(c,1);
    j=cliques(c,2);
    for xi=1:domain_sizes(i)
        for xj=1:domain_sizes(j)
            dl_theta_c(c,xi,xj)=sum_b_c{c}(xi,xj)-lambda*theta_c(c,xi,xj);
        end
    end
end

end