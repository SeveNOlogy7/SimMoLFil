function [ops] = Coupler(name,coupler)
%COUPLER Create two Operation objects in an array for Coupler component
%   name: name of the Coupler Operation
%   coupler: a coupler component
%   A optical coupler usually has 2 input ends and 2 output ends

import simulation.*

ops = [Operation(name,"CouplerOut1",2,...
    struct("coupler",coupler)),...
    Operation(name,"CouplerOut2",2,...
    struct("coupler",coupler))];

end

