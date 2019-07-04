%create_synthetic_scene
%   Randomly orients nb checkerboards of nx by ny points,
%   with the side of the squares of length r.
%   This function returns the scene as a 3xN array of points.
%   There are two optional arguments for the extent of the translation and
%   for turning off rotation of the patterns.
%
%   Usage:
%       points = create_synthetic_scene(6,8,1,5[,0.3[,1]])
%
%   Original code by Simon Donné, January 2017
function [points,base_board,R,t] = create_synthetic_scene(nx,ny,r,nb,translation_fraction, degeneracies)
    if nargin < 5
        translation_fraction = 1.0;
    end
    if nargin < 6
        degeneracies = 0.0;
    end
    
    a = ( (0:nx)-nx/2 )' * ones(1,ny+1);
    b = ones(nx+1,1) * ( (0:ny)-ny/2 );
    base_board = [a(:)';b(:)';zeros(1,numel(a))]*r;
    P = (nx+1)*(ny+1);
    
    points = cell(1,nb);
    R = cell(1,nb);
    t = cell(1,nb);
    for pose = 1:nb
        %we put the checkerboard at twice its height from the camera,
        %give or take a bit
        ti = [0;0;r*(ny+1)*(2-translation_fraction+2*translation_fraction*rand())];
        if rand() > degeneracies
            %random vector for the rotation axis, close to the z-axis
            theta = acos(rand()/4+3/4);
            phi = 2*pi*rand();
            vi = [sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta)];
            %create the rotation matrix using rodrigues' formula
            W =[0 -vi(3) vi(2) ; vi(3) 0 -vi(1) ; -vi(2) vi(1) 0 ];
            psi = pi*rand()-pi/2;
            Ri = eye(3) + sin(psi)*W + 2*sin(psi/2)^2*W^2;
        else
            Ri = eye(3);
        end
        
        %perform the rotation and save the points
        points{pose} = Ri*base_board + repmat(ti,[1,P]);
        R{pose} = Ri;
        t{pose} = ti;
    end
end