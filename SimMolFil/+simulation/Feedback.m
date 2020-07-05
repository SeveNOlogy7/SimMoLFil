function [op] = Feedback()
%FEEDBACK Create a Operation for Feedback
%   Usage:
%       uout = Feedback(u_reccurence,u_in);

import simulation.*

op = Operation("","Feedback",2,...
    struct());

end

