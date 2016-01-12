function [len, p0, p1] = findLS(v)

v = v(:).';
x = diff([0 v 0]);
p0s = [0 find(x == 1)];
p1s = [0 find(x == -1)];
if (length(p0s) ~= length(p1s))
    disp(v);
    fprintf('size of p0s');
    disp(size(p0s));

    fprintf('size of p1s');
    disp(size(p1s));
end
if (length(find(v~=0 & v~=1)) > 0)
    fprintf('XD');
end
[len, idx] = max(p1s - p0s);
p0 = p0s(idx);
p1 = p1s(idx) - 1;
