classdef fractionalDL < handle
    properties (GetAccess=public, SetAccess=protected)
        dl
        freed
    end
    methods % public methods
        function obj = fractionalDL(lag, initValue)
            if nargin == 1
                initValue = 0;
            end
            obj.dl = initFDL(double(lag), double(initValue));
            obj.freed = 0;
        end
        function delete(obj)
            if ~obj.freed
                freeFDL(obj.dl);
                obj.freed = 1;
            end
        end
        function y = process(obj, x)
            if ~obj.freed
                y = procFDL(obj.dl, x);
            else
                y = [];
            end
        end
    end
end