function [op] = Fiber(name,fiber)
%FIBER Create a Operation for Fiber component
%   name: name of the Fiber Operation
%   fiber: a fiber component object, which will eventually be used to
%   calculation pulse propagation
%   Usage:
%       fiber = Fiber(name,fiber);
%       uout = uin + fiber;

import simulation.*

op = Operation(name,"Fiber",1,...
    struct("fiber",fiber));

end

