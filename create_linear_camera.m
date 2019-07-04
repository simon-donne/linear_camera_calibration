%create_linear_camera
%   Creates the intrinsic camera of a linear camera, based on a given focal
%   distance, scan-line height and scanning speed (all in pixels).
%   Generated cameras have a random focal distance within 10% of the given
%   value, have an optical center within 10% of the center of the scan-line
%   and have a scanning speed within 10% of the passed value.
%
%   Usage:
%       K = create_linear_camera(300,160,15)
%
%   Original code by Simon Donné, January 2017
function [K] = create_linear_camera(f,u0,s)
    K = [[f,0,u0];[0,s,0];[0,0,1]];
end