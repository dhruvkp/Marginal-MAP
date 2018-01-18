function [b_i,b_c,s_i,s_c,fval]=qpbo_update_v2(m,gamma,dl_i,dl_c,b_i,b_c,cliques,domain_size,number_variables,number_edges,tolearn,number_tolearn_edges,sample)
%QPBO_UPDATE Summary of this function goes here
%   Detailed explanation goes here

    %%% initializing qpbo inputs %%%
    s_i=b_i;
    s_c=b_c;
    terminalWeights=zeros(length(tolearn),2);
    edgeWeights=zeros(number_tolearn_edges,4);
    edgeList=zeros(2,number_tolearn_edges);
    terminalWeights(:,:)=dl_i{m}(tolearn,:);
    
    cnt=1;
    for c=1:number_edges
        i=cliques(c,1);
        j=cliques(c,2);
        [ti,ii]=ismember(i,tolearn);
        [tj,ij]=ismember(j,tolearn);
        if ti && tj
            edgeList(1,cnt)=ii;
            edgeList(2,cnt)=ij;
            edge=1;
            for xi=1:domain_size
                for xj=1:domain_size
                    edgeWeights(cnt,edge)=dl_c{m,c}(xi,xj);
                    edge=edge+1;
                end
            end
            cnt=cnt+1;
        elseif ti
            terminalWeights(ii,:)=terminalWeights(ii,:)+shiftdim(dl_c{m,c}(:,sample(j)),1);
        elseif tj
            terminalWeights(ij,:)=terminalWeights(ij,:)+shiftdim(dl_c{m,c}(sample(i),:),2);
        end
    end
    edgeList=int32(edgeList);
    terminalWeights=double(terminalWeights);
    edgeWeights=double(edgeWeights);
    
    [nodeBelief, edgeBelief,fval] = QPBO_double_mex(terminalWeights', edgeWeights',edgeList);
    
    for i=1:number_variables
        [ti,ii]=ismember(i,tolearn);
        for xi=1:domain_size
            if ti
                b_i{m}(i,xi)=(1-gamma)*b_i{m}(i,xi)+gamma*nodeBelief(xi,ii);
                s_i{m}(i,xi)=nodeBelief(xi,ii);
            end
        end
    end
    edge=1;
    for c=1:number_edges
        i=cliques(c,1);
        j=cliques(c,2);
        [ti,ii]=ismember(i,tolearn);
        [tj,ij]=ismember(j,tolearn);
        if ti && tj
            cnt=1;
            for xi=1:domain_size
                for xj=1:domain_size
                    b_c{m,c}(xi,xj)=(1-gamma)*b_c{m,c}(xi,xj)+gamma*edgeBelief(cnt,edge);
                    s_c{m,c}(xi,xj)=edgeBelief(cnt,edge);
                    cnt=cnt+1;
                end
            end
            edge=edge+1;
        elseif ti
            for xi=1:domain_size
                for xj=1:domain_size
                    b_c{m,c}(xi,xj)=(1-gamma)*b_c{m,c}(xi,xj)+gamma*nodeBelief(xi,ii)*b_i{m}(j,xj);
                    s_c{m,c}(xi,xj)=nodeBelief(xi,ii)*b_i{m}(j,xj);
                end
            end
        elseif tj
            for xi=1:domain_size
                for xj=1:domain_size
                    b_c{m,c}(xi,xj)=(1-gamma)*b_c{m,c}(xi,xj)+gamma*b_i{m}(i,xi)*nodeBelief(xj,ij);
                    s_c{m,c}(xi,xj)=b_i{m}(i,xi)*nodeBelief(xj,ij);
                end
            end
        end
    end
    
end

