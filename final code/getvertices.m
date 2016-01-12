function [Vq ]=getvertices(gridID, grid)
r = gridID(1);
c = gridID(2);
Vq = zeros(4,2);
Vq (1,:) = grid(r,c,:);
Vq (2,:) = grid(r,c+1,:);
Vq (3,:) = grid(r+1,c,:);
Vq (4,:) = grid(r+1,c+1,:);
end