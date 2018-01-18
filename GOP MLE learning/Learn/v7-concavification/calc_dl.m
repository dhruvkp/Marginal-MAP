function [ dl_i,dl_c,dlp_i,dlp_c ] = calc_dl( b_i,b_c,bp_i,bp_c,w_i,w_i_c,lambda,cliques,order_variables,hidden,output,m )
%CALC_DL calculates delta-L
%   Detailed explanation goes here

dl_i=zeros(size(b_i));
dl_c=zeros(size(b_c));

dlp_i=zeros(size(b_i));
dlp_c=zeros(size(b_c));

[M,number_variables,domain_size]=size(b_i);
number_edges=size(b_c,2);

if m==-1
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
                        if b_c(m,c,xi,xj)~=0
                            dl_c(m,c,xi,xj)=dl_c(m,c,xi,xj) - w_i_c(i,c)- w_i_c(i,c)*log(b_c(m,c,xi,xj));
                        end
                        dlp_c(m,c,xi,xj)=(bp_c(m,c,xi,xj)-b_c(m,c,xi,xj))/lambda  ;
                        if bp_c(m,c,xi,xj)~=0
                            dlp_c(m,c,xi,xj)=dlp_c(m,c,xi,xj)+ w_i_c(i,c)+ w_i_c(i,c)*log(bp_c(m,c,xi,xj));
                        end
                    else
                        dl_c(m,c,xi,xj)=(b_c(m,c,xi,xj)-bp_c(m,c,xi,xj))/lambda ;
                        if b_c(m,c,xi,xj)~=0
                            dl_c(m,c,xi,xj)=dl_c(m,c,xi,xj) - w_i_c(j,c)- w_i_c(j,c)*log(b_c(m,c,xi,xj));
                        end
                        dl_c(m,c,xi,xj)=(bp_c(m,c,xi,xj)-b_c(m,c,xi,xj))/lambda ;
                        if bp_c(m,c,xi,xj)~=0
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
                        dl_i(m,j,xj)=dl_i(m,j,xj)+w_i_c(i,c)-w_i_c(j,c) ;
                        if b_i(m,j,xj)~=0
                            dl_i(m,j,xj)=dl_i(m,j,xj)+(w_i_c(i,c)-w_i_c(j,c))*log(b_i(m,j,xj));
                        end
                    end
                else
                    for xi=1:domain_size
                        dl_i(m,i,xi)=dl_i(m,i,xi)+w_i_c(j,c)-w_i_c(i,c);
                        if b_i(m,i,xi)~=0
                            dl_i(m,i,xi)=dl_i(m,i,xi)+(w_i_c(j,c)-w_i_c(i,c))*log(b_i(m,i,xi));
                        end
                    end
                end
            end
            if ismember(i,horzcat(hidden,output)) && ismember(j,horzcat(hidden,output))
                if find(order_variables==i)<find(order_variables==j) % j is parent of i, j is summed out later than i
                    for xj=1:domain_size
                        dlp_i(m,j,xj)=dlp_i(m,j,xj)-w_i_c(i,c)+w_i_c(j,c) ;
                        if bp_i(m,j,xj)~=0
                            dlp_i(m,j,xj)=dlp_i(m,j,xj)-(w_i_c(i,c)-w_i_c(j,c))*log(bp_i(m,j,xj));
                        end
                    end
                else
                    for xi=1:domain_size
                        dlp_i(m,i,xi)=dlp_i(m,i,xi)-w_i_c(j,c)+w_i_c(i,c);
                        if bp_i(m,i,xi)~=0
                            dlp_i(m,i,xi)=dlp_i(m,i,xi)-(w_i_c(j,c)-w_i_c(i,c))*log(bp_i(m,i,xi));
                        end
                    end
                end
            end
        end
    end
else
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
                    if b_c(m,c,xi,xj)~=0
                        dl_c(m,c,xi,xj)=dl_c(m,c,xi,xj) - w_i_c(i,c)- w_i_c(i,c)*log(b_c(m,c,xi,xj));
                    end
                    dlp_c(m,c,xi,xj)=(bp_c(m,c,xi,xj)-b_c(m,c,xi,xj))/lambda  ;
                    if bp_c(m,c,xi,xj)~=0
                        dlp_c(m,c,xi,xj)=dlp_c(m,c,xi,xj)+ w_i_c(i,c)+ w_i_c(i,c)*log(bp_c(m,c,xi,xj));
                    end
                else
                    dl_c(m,c,xi,xj)=(b_c(m,c,xi,xj)-bp_c(m,c,xi,xj))/lambda ;
                    if b_c(m,c,xi,xj)~=0
                        dl_c(m,c,xi,xj)=dl_c(m,c,xi,xj) - w_i_c(j,c)- w_i_c(j,c)*log(b_c(m,c,xi,xj));
                    end
                    dl_c(m,c,xi,xj)=(bp_c(m,c,xi,xj)-b_c(m,c,xi,xj))/lambda ;
                    if bp_c(m,c,xi,xj)~=0
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
                    dl_i(m,j,xj)=dl_i(m,j,xj)+w_i_c(i,c)-w_i_c(j,c) ;
                    if b_i(m,j,xj)~=0
                        dl_i(m,j,xj)=dl_i(m,j,xj)+(w_i_c(i,c)-w_i_c(j,c))*log(b_i(m,j,xj));
                    end
                end
            else
                for xi=1:domain_size
                    dl_i(m,i,xi)=dl_i(m,i,xi)+w_i_c(j,c)-w_i_c(i,c);
                    if b_i(m,i,xi)~=0
                        dl_i(m,i,xi)=dl_i(m,i,xi)+(w_i_c(j,c)-w_i_c(i,c))*log(b_i(m,i,xi));
                    end
                end
            end
        end
        if ismember(i,horzcat(hidden,output)) && ismember(j,horzcat(hidden,output))
            if find(order_variables==i)<find(order_variables==j) % j is parent of i, j is summed out later than i
                for xj=1:domain_size
                    dlp_i(m,j,xj)=dlp_i(m,j,xj)-w_i_c(i,c)+w_i_c(j,c) ;
                    if bp_i(m,j,xj)~=0
                        dlp_i(m,j,xj)=dlp_i(m,j,xj)-(w_i_c(i,c)-w_i_c(j,c))*log(bp_i(m,j,xj));
                    end
                end
            else
                for xi=1:domain_size
                    dlp_i(m,i,xi)=dlp_i(m,i,xi)-w_i_c(j,c)+w_i_c(i,c);
                    if bp_i(m,i,xi)~=0
                        dlp_i(m,i,xi)=dlp_i(m,i,xi)-(w_i_c(j,c)-w_i_c(i,c))*log(bp_i(m,i,xi));
                    end
                end
            end
        end
    end
end

end