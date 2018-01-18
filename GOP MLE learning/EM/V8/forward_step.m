function alpha = forward_step(model, X)
% Input:
%     model.A: transition matrix (K * k matrix), k no. of hidden states
%     model.B: observation matrix (K * N matrix), N no. of observed states
%     model.S: initial state distribution (1 * k vector)
%     X: ( 1 * N Data vector)
%     
% Output:
%     alpha: (N * K matrix),where every row i holding aplha(zi)

    N = size(X,2);
    K = size(model.S,2);
    alpha = zeros(N, K);
    % intialize alpha(z1)
    alpha(1,:) = exp(model.S) .* exp(model.B(:, X(1)))';

    % Calculate the values of the remaining alphas
    for i = 2:N
        alpha(i,:) = exp(model.B(:, X(i)))' .* (alpha(i-1, :) * exp(model.A) );
    end
end