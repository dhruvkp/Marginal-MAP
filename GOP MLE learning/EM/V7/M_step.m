function [new_model, llh] = M_step(old_model, X)
% Input:
%     old_model.A: transition matrix (K * k matrix), k no. of hidden states
%     old_model.B: observation matrix (K * N matrix), N no. of observed states
%     old_model.S: initial state distribution (1 * k vector)
%     X: ( 1 * N Data vector)
%     
% Output:
%    new_model have same structure as old_model
%    llh: log likelihood

% intialize our variables
N = size(X,2);
K = size(old_model.S,2);
% call E_step to calculate gamma and xi
[gamma, xi, llh] = E_step(old_model, X);

% set new_model.S
new_model.S = gamma(1,:)./ sum(gamma(1,:));

% set new_model.A
% intialize it
new_model.A = xi{1,1};
% for every entry in xi add it to A
for i =2:N-1
new_model.A = new_model.A + xi{i,1};
end

% replace all zero vectors with ones
zero_entries = all(new_model.A==0,2);
new_model.A(zero_entries,:) = repmat(ones(1,K),nnz(zero_entries) ,1);

% normalize A
new_model.A = new_model.A ./ sum(new_model.A,2);

%set new_model.B
repeated_data = repmat(X, size(old_model.B,2),1);
occurrence =  cell2mat(arrayfun(@(i) all(repeated_data(i,:)==i,1),1:size(repeated_data,1),'un',0)');
new_model.B = (occurrence * gamma)';
new_model.B = new_model.B ./ sum(new_model.B, 2);
end