%linear_camera_projection
%   Projects a 3xK array of points using a linear camera matrix
%   The points are assumed to lie in the camera space
%   Usage:
%       uv = linear_camera_projection(K,XYZ)
%
%   Original code by Simon Donné, January 2017
function [projected] = linear_camera_projection(camera,points)
    projected = [(camera(1,:)*points)./(camera(3,:)*points); ...
                 camera(2,:)*points];
end