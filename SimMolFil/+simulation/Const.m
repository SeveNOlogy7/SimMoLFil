function [m] = Const(value)
%COUPLER Create a model for Const signal
%   Const Models need no name
%   value: value hold by the model

import simulation.*

m = Model(Operation("", "Const", 0,...
    struct("value",value)),[]);

end

