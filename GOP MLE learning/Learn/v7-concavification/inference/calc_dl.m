function [ dl_i,dl_c ] = calc_dl( b_i,b_c,theta_i,theta_c,w_i,w_i_c,cliques,order_variables,hidden,output,domain_sizes )
%CALC_DL calculates delta-L
%   Detailed explanation goes here

dl_i=b_i;
dl_c=b_c;

number_variables=size(b_i,1);
number_edges=size(b_c,1);
for i=1:number_variables
    for xi=1:domain_sizes(i)
        dl_i(i,xi)=-theta_i(i,xi);
        if b_i(i,xi)~=0 && (ismember(i,hidden) || ismember(i,output))
            dl_i(i,xi)=dl_i(i,xi)+ w_i(i)+ w_i(i)*log(b_i(i,xi));
        end
    end
end

for c=1:number_edges
    i=cliques(c,1);
    j=cliques(c,2);
    for xi=1:domain_sizes(i)
        for xj=1:domain_sizes(j)
            dl_c(c,xi,xj)=-theta_c(i,j,xi,xj);
            if find(order_variables==i)<find(order_variables==j) % j is parent of i
                if b_c(c,xi,xj)~=0
                    dl_c(c,xi,xj)=dl_c(c,xi,xj)+ w_i_c(i,c)+ w_i_c(i,c)*log(b_c(c,xi,xj));
                end
            else
                if b_c(c,xi,xj)~=0
                    dl_c(c,xi,xj)=dl_c(c,xi,xj) + w_i_c(j,c)+ w_i_c(j,c)*log(b_c(c,xi,xj));
                end
            end
        end
    end
end

for c=1:number_edges
    i=cliques(c,1);
    j=cliques(c,2);
    if ismember(i,horzcat(hidden,output)) && ismember(j,horzcat(hidden,output))
        if find(order_variables==i)<find(order_variables==j) % j is parent of i, j is summed out later than i
            for xj=1:domain_sizes(j)
                dl_i(j,xj)=dl_i(j,xj)-w_i_c(i,c)+w_i_c(j,c) ;
                if b_i(j,xj)~=0
                    dl_i(j,xj)=dl_i(j,xj)-(w_i_c(i,c)-w_i_c(j,c))*log(b_i(j,xj));
                end
            end
        else
            for xi=1:domain_sizes(i)
                dl_i(i,xi)=dl_i(i,xi)-w_i_c(j,c)+w_i_c(i,c);
                if b_i(i,xi)~=0
                    dl_i(i,xi)=dl_i(i,xi)-(w_i_c(j,c)-w_i_c(i,c))*log(b_i(i,xi));
                end
            end
        end
    end
end

end