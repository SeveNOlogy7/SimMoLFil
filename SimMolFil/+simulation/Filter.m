function [op] = Filter(name, filter)
%FILTER Summary of this function goes here
%   Detailed explanation goes here

import simulation.*

op = Operation(name,"Filter",1,...
    struct("filter",filter));

end

