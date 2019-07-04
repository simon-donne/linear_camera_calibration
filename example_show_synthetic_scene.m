%show_synthetic_example
%   Visualizes the generated synthetic scenes
%   Usage:
%       show_synthetic_example()
%
%   Original code by Simon Donné, January 2017
function [] = example_show_synthetic_scene()


close all

nx = 6;
ny = 8;
r = 1;
nb = 5;

rng(1)
points = create_synthetic_scene(nx,ny,r,nb,0.5,0);

bb = [ [Inf,-Inf];[Inf,-Inf];[Inf,-Inf]];

figure, axes, hold all, axis equal, axis vis3d
for b = 1:nb
    scatter3(points{b}(1,:),points{b}(2,:),points{b}(3,:))
    bbb = [min(points{b},[],2),max(points{b},[],2)];
    bb(:,1) = min(bb(:,1),bbb(:,1));
    bb(:,2) = max(bb(:,2),bbb(:,2));
end
plot_box(bb);

end

function [] = plot_box(bb)
    color = 'k';
    
    plot3([bb(1,1),bb(1,2)],[bb(2,1),bb(2,1)],[bb(3,1),bb(3,1)],color)
    plot3([bb(1,1),bb(1,1)],[bb(2,1),bb(2,2)],[bb(3,1),bb(3,1)],color)
    plot3([bb(1,1),bb(1,1)],[bb(2,1),bb(2,1)],[bb(3,1),bb(3,2)],color)
    
    plot3([bb(1,2),bb(1,1)],[bb(2,2),bb(2,2)],[bb(3,1),bb(3,1)],color)
    plot3([bb(1,2),bb(1,2)],[bb(2,2),bb(2,1)],[bb(3,1),bb(3,1)],color)
    plot3([bb(1,2),bb(1,2)],[bb(2,2),bb(2,2)],[bb(3,1),bb(3,2)],color)
    
    plot3([bb(1,1),bb(1,2)],[bb(2,2),bb(2,2)],[bb(3,2),bb(3,2)],color)
    plot3([bb(1,1),bb(1,1)],[bb(2,2),bb(2,1)],[bb(3,2),bb(3,2)],color)
    plot3([bb(1,1),bb(1,1)],[bb(2,2),bb(2,2)],[bb(3,2),bb(3,1)],color)
    
    plot3([bb(1,2),bb(1,1)],[bb(2,1),bb(2,1)],[bb(3,2),bb(3,2)],color)
    plot3([bb(1,2),bb(1,2)],[bb(2,1),bb(2,2)],[bb(3,2),bb(3,2)],color)
    plot3([bb(1,2),bb(1,2)],[bb(2,1),bb(2,1)],[bb(3,2),bb(3,1)],color)
end