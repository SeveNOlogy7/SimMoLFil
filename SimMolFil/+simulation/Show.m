function [op] = Show(options)
%SHOW create a show operation to display a model output
%   options: struct that contain show parameters
%       After creation, Show operator s can be called by m = s(m)
%       Show operator should return same model as it input
arguments
    options struct
end

import simulation.*

op = Operation("", "Show", 1,...
    options);

end

