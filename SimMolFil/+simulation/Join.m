function [op] = Join()
%JOIN Create a Operation for Feedback
%   Usage:
%       u_main = Join(u_main+u_side);
%   Extende to u_main = Join(u_main+u_side1+...) in the future

import simulation.*

op = Operation("","Join",2,...
    struct());

end

