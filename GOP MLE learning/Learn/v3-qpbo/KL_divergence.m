function [ KLD ] = KL_divergence( A,output,input,hidden,theta_true,theta_learned)
%KL_DIVERGENCE Summary of this function goes here
%   Detailed explanation goes here
number_variables=size(A,1);
KLD=0;

for x=binary_combinations(length(input))'
    px1=0;
    px2=0;
    for y=binary_combinations(length(output))'
        for h=binary_combinations(length(hidden))'
            v=zeros(1,number_variables);
            v(output)=y;
            v(input)=x;
            v(hidden)=h;
            px1=px1+calc_p_all_x(A,v,theta_true);
            px2=px2+calc_p_all_x(A,v,theta_learned);
           
        end
    end
    for y=binary_combinations(length(output))'
        p1=0;
        p2=0;
        for h=binary_combinations(length(hidden))'
            v=zeros(1,number_variables);
            v(output)=y;
            v(input)=x;
            v(hidden)=h;
            p1=p1+calc_p_all_x(A,v,theta_true);
            p2=p2+calc_p_all_x(A,v,theta_learned);
        end
        p=p1/px1;
        q=p2/px2;
        KLD=KLD+p*log(p/q);
    end
end

end

