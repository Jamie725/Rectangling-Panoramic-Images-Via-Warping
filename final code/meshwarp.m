function [ outputimg ] = meshwarp(imgin,Vlocal, Vglobal)

%% perform warping based on Vwarp and Vopt 
[rows,cols,color] = size(imgin);
outputimg = zeros(rows,cols,3);
outputimgcount = zeros(rows,cols,3);
Vwarp = Vlocal;
Vopt = Vglobal;
[ygridNum , xgridNum, xy] = size(Vwarp);
quadrows = ygridNum -1;
quadcols = xgridNum -1;
for i = 1:quadrows
    for j = 1:quadcols
        Vq=getvertices([i j], Vopt);
        Vq = fliplr(Vq);
        Vq = reshape(Vq',[],1); % Vq should be x y
        Vo=getvertices([i j], Vwarp);
        Vo = fliplr(Vo);
        Vo = reshape(Vo',[],1); % Vq should be x y
        xlength = max(Vq(1:2:8)) -min(Vq(1:2:8));
        ylength = max(Vq(2:2:8)) -min(Vq(2:2:8));
        t2nstep = 1/(xlength*4);
        t1nstep = 1/(ylength*4);
        for t1n = 0:t1nstep:1
            for t2n = 0:t2nstep:1
                v1w=1-t1n-t2n+t1n*t2n;
                v2w=t2n-t1n*t2n;
                v3w=t1n-t1n*t2n;
                v4w=t1n*t2n;
                T = [ v1w 0 v2w 0 v3w 0 v4w 0; ...
                        0 v1w 0 v2w 0 v3w 0 v4w ];
                pout = round(T*Vq);
                pref = round(T*Vo); 
                try  % Vq may contain negetive index
                outputimg(pout(2),pout(1),:) = outputimg(pout(2),pout(1),:)+ imgin(pref(2),pref(1),:);
                outputimgcount (pout(2),pout(1),:) = outputimgcount (pout(2),pout(1),:) +1;
                catch
                    
                end
            end
        end
        figure(4); imshow(outputimg./outputimgcount);
    end
end
outputimg= outputimg./outputimgcount;

% imwrite(outputimg,[datapath 'global.png']);
% 
% gridmask = drawGridmask(Vopt,  rows,cols);
% imageGrided2 = drawGrid( gridmask, outputimg);
% figure(5); imshow(imageGrided2);
% imwrite(imageGrided2,[datapath 'global_mesh.png']);


end

