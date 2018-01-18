function [samples]=gibbs_sampler_mrf_with_edge_parameters(A,theta_c,domain_sizes,burnin,number_samples)
% input: adjacency matrix, clique parameters, burn-in iterations, number of
% samples required
% output: samples
number_variables=size(A,1);
order=1:number_variables;
x=ones(1,number_variables);
while calc_p_all_x(A,x,theta_c)==0
    for i=1:number_variables
        x(i)=randi(domain_sizes(i),1);
    end
end
samples=zeros(number_samples,number_variables);
for cnt=1:(burnin+number_samples)
    for j=order
        p=zeros(domain_sizes(j),1);
        x_domain=1:domain_sizes(j);
        for x_j=x_domain
            value=1;
            for i=1:number_variables
                for k=1:number_variables
                    if A(i,k)>0 && i>k
                        if i==j
                            value=value*exp(theta_c(i,k,x_j,x(k)));
                        elseif k==j
                            value=value*exp(theta_c(i,k,x(i),x_j));
                        else
                            value=value*exp(theta_c(i,k,x(i),x(k)));
                        end
                    end
                end
            end
            p(x_j)=value;
        end
        p=p/norm(p,1);
        p=cumsum(p);
        r=rand();
        for i=x_domain
            if p(i)>r
                x(j)=i;
                break;
            end
        end
    end
    if cnt>burnin
        samples(cnt-burnin,:)=x;
    end
end