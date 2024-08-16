function output = anglediff(angle1,angle2)

% calculate the angle difference betwen angle 1 and angle 2
% please make sure all the inputs are in the range [0 360];
% By Gao Yuan 2023-11-01 19:32
    normDeg = mod(angle1-angle2-180,360)+180;
    output = min(360-normDeg, normDeg);
 end
 
