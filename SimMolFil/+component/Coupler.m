classdef Coupler
    %COUPLER simple four-end optical coupler
    %   u1i       u2i
    %      |     |
    %        |-|
    %        |-|       
    %      |     |
    %   u1o       u2o
    
    properties
        rho     % couple ratio (u1i,u2i)->(u1o,u2o), u1o = sqrt(rho)*u1i
    end
    
    methods
        function obj = Coupler(rho)
            %COUPLER Construct an instance of this class
            arguments
                rho {mustBeGreaterThanOrEqual(rho,0),mustBeLessThanOrEqual(rho,1)}
            end
            obj.rho = rho;
        end
        
        function u1o = couplerOut1(obj,u1i,u2i)
            u1o = sqrt(obj.rho)*u1i + 1i*sqrt(1-obj.rho)*u2i;
        end
        
        function u2o = couplerOut2(obj,u1i,u2i)
            u2o = 1i*sqrt(1-obj.rho)*u1i + sqrt(obj.rho)*u2i;
        end
    end
end

