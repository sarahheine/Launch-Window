% A function is created which reduces angles larger than 2*Pi such that
% they are smaller then 2*Pi.

% Micro-X Launch Window Calculations
% Function Reduce
% Version 1.0
    
function [SmallAngle] = Reduce(LargeAngle)
    while (LargeAngle > 2*pi)
        LargeAngle = LargeAngle - 2*pi;
    end
    SmallAngle = LargeAngle;
end
