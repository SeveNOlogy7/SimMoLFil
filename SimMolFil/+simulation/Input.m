function [m] = Input(name)
%INPUT Create a Model object that asks for input
%   name: name of the input operation
%   some input will eventually be put into a variable called name

import simulation.*

m = Model(Operation(name, "Input", 0, struct()),[]);

end

