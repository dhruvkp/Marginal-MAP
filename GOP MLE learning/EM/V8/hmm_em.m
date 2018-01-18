function model = hmm_em(hidden, observed,training_data, parameters)
% Input:
%     hidden: number of hidden values
%     hidden: number of observed values
%     training_data: cell array of size N * 1, where N is the number of
%                    training instances.
%     parameters: set of parameters (max iteration number, tolerance, LR)
%    
%     
% Output:
%     model.A: transition matrix (K * k matrix), k no. of hidden states
%     model.B: observation matrix (K * N matrix), N no. of observed states
%     model.S: initial state distribution (1 * k vector)

model = initialize(hidden, observed);
i = 1;
llh_old = 1;
llh_new = 0;
while i<parameters.iterations %&& abs(llh_old - llh_new) > parameters.tolerance
    i = i+1;
    llh_old = llh_new;
    temp_model = M_step(model, training_data{1,1});
    for j = 2:size(training_data,1)
        [new_model, llh] = M_step(model, training_data{j,1});
        llh_new = llh_new + log(llh);
        temp_model.A = temp_model.A + new_model.A;
        temp_model.B = temp_model.B + new_model.B;
        temp_model.S = temp_model.S + new_model.S;
    end
    model.A = model.A +  parameters.LR* ((temp_model.A / size(training_data,1))- (parameters.lambda* sum(sum(model.A)) ));
    model.B = model.B + parameters.LR* ((temp_model.B / size(training_data,1))- (parameters.lambda* sum(sum(model.B)) ));
    model.S = model.S + parameters.LR* ((temp_model.S / size(training_data,1))- (parameters.lambda* sum(model.S) ));
    llh_new = llh_new / size(training_data,1);
%     model.A = model.A - mean(mean(model.A));
%     model.B = model.B - mean(mean(model.B));
%     model.S = model.S - mean(mean(model.S));
end
end