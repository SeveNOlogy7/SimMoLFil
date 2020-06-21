function [op] = Switch(name,condition)
%SWITCH Create a operator for switching signals (models)
%   name: name of the switch model
%   condition: string that hold the condtion. 
%       After creation, Switch operator s can be called by s(m1,m2)
%       In implementation, if condition is true, m = m1
%   
import simulation.*

op = Operation(name, "Switch", 2,...
    struct("condition",condition));

end

