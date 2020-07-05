function [op] = Filter(name, filter)
%FILTER Create a Filter Operation
%   Usage:
%       filter = Filter(name, filter);
%       uout = uin + filter;

import simulation.*

op = Operation(name,"Filter",1,...
    struct("filter",filter));

end

