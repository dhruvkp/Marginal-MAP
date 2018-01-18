function [gamma, xi,llh] = E_step(model, X)
% Input:
%     model.A: transition matrix (K * k matrix), k no. of hidden states
%     model.B: observation matrix (K * N matrix), N no. of observed states
%     model.S: initial state distribution (1 * k vector)
%     X: ( 1 * N Data vector)
%     
% Output:
%     gamma: (N * K matrix), holding the marginal posterior distribution of the latent variables   
%     xi: a structure of N-1 (K * K matrices), each holds the joint posterior distribution of two
%         successive latent variables
%    llh: log likelihood

% intialize our variables
N = size(X,2);
K = size(model.S,2);
gamma = zeros(N, K);
xi = cell(N-1,1);

% call forward_step to calculate alpha
alpha = forward_step(model, X);

% call backward step to calculate beta
beta = backward_step(model, X);

% calculate P(X)
%Pro_of_X = exp(sum(alpha(N,:)));
llh = sum(alpha(N,:));
% calculate gamma
gamma = (alpha .* beta)/llh;%./Pro_of_X;

% calculate xi

for i = 1:N-1
     xi{i,1}= (alpha(i,:)' .* exp(model.B(: , X(i+1)))' .* exp(model.A(: , :)) .* beta(i+1 , :))/llh;%/ Pro_of_X;
end
 

end