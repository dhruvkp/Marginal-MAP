function beta = backward_step(model, X)
% Input:
%     model.A: transition matrix (K * k matrix), k no. of hidden states
%     model.B: observation matrix (K * N matrix), N no. of observed states
%     model.S: initial state distribution (1 * k vector)
%     X: ( 1 * N Data vector)
%     
% Output:
%     beta: (N * K matrix),where every row i holding beta(zi)

    N = size(X,2);
    K = size(model.S,2);
    beta = ones(N, K);
    % intialize beta(z1)
    beta(1,:) = ones(1,K);

    % Calculate the values of the remaining betas
    for i = N-1:-1:1
        beta(i,:) = (exp(model.B(:, X(i+1)))' .* beta(i+1, :)) * exp(model.A)';
    end
end