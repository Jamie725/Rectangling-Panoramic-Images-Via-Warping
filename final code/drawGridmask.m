function [ gridmask ] = drawGridmask( V, rows,cols )
ygrid = round(V(:,:,1));
xgrid = round(V(:,:,2));
xgridN=size(V,2);
ygridN=size(V,1);
gridmask = false(rows,cols);
for x=1:xgridN
    for y=1:ygridN
        try
            selfpoint = [ ygrid(y,x) xgrid(y,x)];
            uppoint = [ ygrid(y-1,x) xgrid(y-1,x)];
            for i=uppoint(1):selfpoint(1)
                m=(selfpoint(2)-uppoint(2))/(selfpoint(1)-uppoint(1));
                gridmask(i, uppoint(2) + round(m*(i-uppoint(1))) )=1;
            end
        catch
        end
        try
            selfpoint = [ ygrid(y,x) xgrid(y,x)];
            leftpoint = [ ygrid(y,x-1) xgrid(y,x-1)];
            for j=leftpoint(2):selfpoint(2)
                m=(selfpoint(1)-leftpoint(1))/(selfpoint(2)-leftpoint(2));
                gridmask( leftpoint(1) + round(m*(j-leftpoint(2))), j)=1;
            end
        catch
        end
    end
end

end

