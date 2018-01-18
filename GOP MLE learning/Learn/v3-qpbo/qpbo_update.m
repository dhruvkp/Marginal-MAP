function [b_i,b_c,s_i,s_c]=qpbo_update(m,gamma,dl_i,dl_c,b_i,b_c,cliques,domain_size,number_variables,number_edges,tolearn,number_tolearn_edges)
%QPBO_UPDATE Summary of this function goes here
%   Detailed explanation goes here

    %%% initializing qpbo inputs %%%
    s_i=b_i;
    s_c=b_c;
    terminalWeights=zeros(length(tolearn),2);
    edgeWeights=zeros(number_tolearn_edges,6);
    terminalWeights(:,:)=dl_i(m,tolearn,:);
    
    cnt=1;
    for c=1:number_edges
        i=cliques(c,1);
        j=cliques(c,2);
        [ti,ii]=ismember(i,tolearn);
        [tj,ij]=ismember(j,tolearn);
        if ti && tj
            edgeWeights(cnt,1)=ii;
            edgeWeights(cnt,2)=ij;
            edge=3;
            for xi=1:domain_size
                for xj=1:domain_size
                    edgeWeights(cnt,edge)=dl_c(m,c,xi,xj);
                    edge=edge+1;
                end
            end
            cnt=cnt+1;
        end
    end
    edgeWeights
   
    [lowerBound, labels] = qpboMex(terminalWeights, edgeWeights);
    labels
    
    for i=1:number_variables
        [ti,ii]=ismember(i,tolearn);
        for xi=1:domain_size
            if ti
                if labels(ii)<0
                    b_i(m,i,xi)=(1-gamma)*b_i(m,i,xi)+gamma*1/2;
                    s_i(m,i,xi)=1/2;
                elseif labels(ii)==(xi-1)
                    b_i(m,i,xi)=(1-gamma)*b_i(m,i,xi)+gamma;
                    s_i(m,i,xi)=1;
                else
                    b_i(m,i,xi)=(1-gamma)*b_i(m,i,xi);
                    s_i(m,i,xi)=0;
                end
            end
        end
    end
    cnt=1;
    for c=1:number_edges
        i=cliques(c,1);
        j=cliques(c,2);
        [ti,ii]=ismember(i,tolearn);
        [tj,ij]=ismember(j,tolearn);
        if ti && tj
            if labels(ii)<0 && labels(ij)<0
                if edgeWeights(cnt,3)+edgeWeights(cnt,6) > edgeWeights(cnt,4)+edgeWeights(cnt,5)
                    b_c(m,c,0,1)=0.5;
                    b_c(m,c,1,0)=0.5;
                else
                    b_c(m,c,0,0)=0.5;
                    b_c(m,c,1,1)=0.5;
                end
            else
                for xi=1:domain_size
                    for xj=1:domain_size
                        b_c(m,c,xi,xj)=(1-gamma)*b_c(m,c,xi,xj)+gamma*(((xi-1)==labels(ii))+(1/2)*(labels(ii)<0))*(((xj-1)==labels(ij))+(1/2)*(labels(ij)<0));
                        s_c(m,c,xi,xj)=(((xi-1)==labels(ii))+(1/2)*(labels(ii)<0))*(((xj-1)==labels(ij))+(1/2)*(labels(ij)<0));
                    end
                end
            end
            cnt=cnt+1;
        elseif ti
            for xi=1:domain_size
                for xj=1:domain_size
                    b_c(m,c,xi,xj)=(1-gamma)*b_c(m,c,xi,xj)+gamma*(((xi-1)==labels(ii))+(1/2)*(labels(ii)<0))*b_i(m,j,xj);
                    s_c(m,c,xi,xj)=(((xi-1)==labels(ii))+(1/2)*(labels(ii)<0))*b_i(m,j,xj);
                end
            end
        elseif tj
            for xi=1:domain_size
                for xj=1:domain_size
                    b_c(m,c,xi,xj)=(1-gamma)*b_c(m,c,xi,xj)+gamma*b_i(m,i,xi)*(((xj-1)==labels(ij))+(1/2)*(labels(ij)<0));
                    s_c(m,c,xi,xj)=b_i(m,i,xi)*(((xj-1)==labels(ij))+(1/2)*(labels(ij)<0));
                end
            end
        end
    end
    
end

