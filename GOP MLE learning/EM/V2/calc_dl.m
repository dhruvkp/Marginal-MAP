function [ dl_i,dl_c,dl_pa_i,dlp_i,dlp_c,dlp_pa_i ] = calc_dl( b_i,b_c,b_pa_i,bp_i,bp_c,bp_pa_i,w_i,w_i_c,lambda,cliques,order_variables,hidden,output )
%CALC_DL calculates delta-L
%   Detailed explanation goes here

dl_i=zeros(size(b_i));
dl_c=zeros(size(b_c));
dl_pa_i=zeros(size(b_pa_i));

dlp_i=zeros(size(b_i));
dlp_c=zeros(size(b_c));
dlp_pa_i=zeros(size(b_pa_i));

[M,number_variables,domain_size]=size(b_i);
number_edges=size(b_c,2);

for m=1:M
    for i=1:number_variables
        for xi=1:domain_size
            dl_i(m,i,xi)=(b_i(m,i,xi)-bp_i(m,i,xi))/lambda;
            if b_i(m,i,xi)~=0 && ismember(i,hidden)
                dl_i(m,i,xi)=dl_i(m,i,xi) - w_i(i)- w_i(i)*log(b_i(m,i,xi));
            end
            dlp_i(m,i,xi)=(bp_i(m,i,xi)-b_i(m,i,xi))/lambda  ;
            if bp_i(m,i,xi)~=0 && (ismember(i,hidden) || ismember(i,output))
                dlp_i(m,i,xi)=dlp_i(m,i,xi)+ w_i(i)+ w_i(i)*log(bp_i(m,i,xi));
            end
        end
    end
    
    for c=1:number_edges
        i=cliques(c,1);
        j=cliques(c,2);
        for xi=1:domain_size
            for xj=1:domain_size
                if find(order_variables==i)<find(order_variables==j) % j is parent of i
                    dl_c(m,c,xi,xj)=(b_c(m,c,xi,xj)-bp_c(m,c,xi,xj))/lambda ;
                    if b_c(m,c,xi,xj)~=0 &&  ismember(i,hidden)
                        dl_c(m,c,xi,xj)=dl_c(m,c,xi,xj) - w_i_c(i,c)- w_i_c(i,c)*log(b_c(m,c,xi,xj));
                    end
                    dlp_c(m,c,xi,xj)=(bp_c(m,c,xi,xj)-b_c(m,c,xi,xj))/lambda  ;
                    if bp_c(m,c,xi,xj)~=0 && (ismember(i,hidden) || ismember(i,output))
                        dlp_c(m,c,xi,xj)=dlp_c(m,c,xi,xj)+ w_i_c(i,c)+ w_i_c(i,c)*log(bp_c(m,c,xi,xj));
                    end
                else
                    dl_c(m,c,xi,xj)=(b_c(m,c,xi,xj)-bp_c(m,c,xi,xj))/lambda ;
                    if b_c(m,c,xi,xj)~=0 &&  ismember(j,hidden)
                        dl_c(m,c,xi,xj)=dl_c(m,c,xi,xj) - w_i_c(j,c)- w_i_c(j,c)*log(b_c(m,c,xi,xj));
                    end
                    dl_c(m,c,xi,xj)=(bp_c(m,c,xi,xj)-b_c(m,c,xi,xj))/lambda ;
                    if bp_c(m,c,xi,xj)~=0 && (ismember(j,hidden) || ismember(j,output))
                        dlp_c(m,c,xi,xj)=dlp_c(m,c,xi,xj) + w_i_c(j,c)+ w_i_c(j,c)*log(bp_c(m,c,xi,xj));
                    end
                end
            end
        end
    end
    
    for c=1:number_edges
        i=cliques(c,1);
        j=cliques(c,2);
        if ismember(i,hidden) && ismember(j,hidden)
            if find(order_variables==i)<find(order_variables==j) % j is parent of i, j is summed out later than i
                for xj=1:domain_size
                    dl_pa_i(m,c,xj)=w_i_c(i,c)-w_i_c(j,c) ;
                    if b_pa_i(m,c,xj)~=0
                        dl_pa_i(m,c,xj)=dl_pa_i(m,c,xj)+(w_i_c(i,c)-w_i_c(j,c))*log(b_pa_i(m,c,xj));
                    end
                end
            else
                for xi=1:domain_size
                    dl_pa_i(m,c,xi)=w_i_c(j,c)-w_i_c(i,c) +(w_i_c(j,c)-w_i_c(i,c))*log(b_pa_i(m,c,xi));
                    if b_pa_i(m,c,xi)~=0
                        dl_pa_i(m,c,xi)=dl_pa_i(m,c,xi)+(w_i_c(j,c)-w_i_c(i,c))*log(b_pa_i(m,c,xi));
                    end
                end
            end
        end
        if ismember(i,horzcat(hidden,output)) && ismember(j,horzcat(hidden,output))
            if find(order_variables==i)<find(order_variables==j) % j is parent of i, j is summed out later than i
                for xj=1:domain_size
                    dlp_pa_i(m,c,xj)=-w_i_c(i,c)+w_i_c(j,c) ;
                    if bp_pa_i(m,c,xj)~=0
                        dlp_pa_i(m,c,xj)=dlp_pa_i(m,c,xj)-(w_i_c(i,c)-w_i_c(j,c))*log(bp_pa_i(m,c,xj));
                    end
                end
            else
                for xi=1:domain_size
                    dlp_pa_i(m,c,xi)=-w_i_c(j,c)+w_i_c(i,c) -(w_i_c(j,c)-w_i_c(i,c))*log(bp_pa_i(m,c,xi));
                    if bp_pa_i(m,c,xi)~=0
                        dlp_pa_i(m,c,xi)=dlp_pa_i(m,c,xi)-(w_i_c(j,c)-w_i_c(i,c))*log(bp_pa_i(m,c,xi));
                    end
                end
            end
        end
    end
end


end