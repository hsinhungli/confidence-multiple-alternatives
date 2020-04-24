%Experiment to fit (1,2 or 3)
exp = 1;

%Subject to fit (There are 13 subjects in Experiment 1, 11 subjects in Experiment 2, and 11 subjects in Experiment 3)
subject = 1;

%Model to fit
%1: Max model (with both sensory nosie and inference noise)
%2: Diff model
%3: Ent model
model = 1;

output = fitmodel_on_individual(exp, subject, model);