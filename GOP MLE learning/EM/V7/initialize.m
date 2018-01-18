function model = initialize(hidden, observed)
model.S = rand(1, hidden);
model.S = model.S / sum(model.S);
model.A = rand(hidden, hidden);
model.A = model.A ./ sum(model.A,2);
model.B = rand(hidden, observed);
model.B = model.B ./ sum(model.B, 2);
end