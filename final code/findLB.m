function [direction, p0, p1] = findLB(mask)

bnds = cell(4, 1);

bnds{1} = mask(1, :);    % up
bnds{2} = mask(end, :);  % down
bnds{3} = mask(:, 1);    % left
bnds{4} = mask(:, end);  % right

direction = 0;
longest = 0;
p0 = 0;
p1 = 0;
for i = 1:4
	[len, q0, q1] = findLS(bnds{i});
	if (len > longest)
        longest = len;
		direction = i;
		p0 = q0;
		p1 = q1;
	end
end
