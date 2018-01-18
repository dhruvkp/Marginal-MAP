function [ KLD ] = KL_divergence( A,output,input,hidden,theta_true,theta_learned)
%KL_DIVERGENCE Summary of this function goes here
%   Detailed explanation goes here
number_variables=size(A,1);
number_edges=sum(sum(A))/2;
theta_max=mean(abs(theta_learned(:)));
theta_learned=theta_learned-theta_max;

KLD=0;

for x=binary_combinations(length(input))'
    px=0;
    qx=0;
    for y=binary_combinations(length(output))'
        for h=binary_combinations(length(hidden))'
            v=zeros(1,number_variables);
            v(output)=y;
            v(input)=x;
            v(hidden)=h;
            px=px+calc_p_all_x(A,v,theta_true);
            qx=qx+calc_p_all_x(A,v,theta_learned);
           
        end
    end
    for y=binary_combinations(length(output))'
        pxy=0;
        qxy=0;
        for h=binary_combinations(length(hidden))'
            v=zeros(1,number_variables);
            v(output)=y;
            v(input)=x;
            v(hidden)=h;
            pxy=pxy+calc_p_all_x(A,v,theta_true);
            qxy=qxy+calc_p_all_x(A,v,theta_learned);
        end
        p=pxy/px;
        q=qxy/qx;
        KLD=KLD+p*log(p/q);
    end
end

end

