function [ lambda,KLD ] = parameter_tuning( A,x_domain,input,output,hidden,samples,theta_c,number_runs,solver )
%PARAMTER_TUNING Summary of this function goes here
%   Detailed explanation goes here

lambda_range=[5*10^-6,5*10^-5,5*10^-4,0.005,0.5,1];
if ~strcmp(solver,'EM')
    lambda_range=[1,10,100,1000,10000];
end


KLD_values=zeros(number_runs,length(lambda_range));


cnt=1;
for lambda=lambda_range   
    for n=1:number_runs
        if strcmp(solver,'EM')
            theta_c_learned = EM( A,size(theta_c),hidden,output,input,samples,lambda);
        else
            theta_c_learned=learn(A,input,output,hidden,x_domain,samples,lambda,solver);
        end
        KLD_values(n,cnt)=KL_divergence( A,output,input,hidden,theta_c,theta_c_learned);
        [n;lambda]
    end
    cnt=cnt+1;
end
% 
% figure
% % errorbar(sample_sizes,mean(KLD_values,1),min(KLD_values,[],1)-mean(KLD_values,1),max(KLD_values,[],1)-mean(KLD_values,1),'-g.')
% % hold on
% errorbar(log(lambda_range),mean(KLD_values,1),min(KLD_values,[],1)-mean(KLD_values,1),max(KLD_values,[],1)-mean(KLD_values,1),'-b.')
% ylabel('KL divergence of p(Y|X,theta)');
% lambda_symbol = char(955);     % lambda unicode
% xlabel(['Log(' lambda_symbol ')']);
% legend(solver);
KLD_avg=mean(KLD_values,1);
[~,idx]=min(KLD_avg);
lambda=lambda_range(idx);
KLD=KLD_values(:,idx);
end

