function [m] = Recurrence(name, initial_value, condition)
%VARIABLE Create a Model presents a Recuurence Point
%   name: name of the Recurence Point in this model
%   initial_value: the initial value of the recurrence node (Model/const)
%   condition: if condition is satified, feedback to Recurrence
%       Once created, Recuurence v can be used with feedback operatior f
%       and another variable s. By writting v = f(v, s), the programe will
%       go to the statment when v is created and update the value of the
%       Recurrence point v until the condition turn false

import simulation.*

m = Model(Operation(name, "Recurrence", 1,...
    struct("condition",condition)),Model(initial_value));
end

