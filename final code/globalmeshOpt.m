function [ Vlocal, Vglobal ] = globalmeshOpt( origImg,mask,dispMap,name )
%GLOBALWARPING Summary of this function goes here
%   Detailed explanation goes here
%datapath = name;

%load(['data\' datapath]);
%load('mask.mat');
% Need mask and origImg

imgin =origImg;

%[ dispMap, outImg ] = localwarping( imgin, mask );
%load('dispMap.mat'); % my dispmap implementation
imagepathforlsd='imglsd.png';
imwrite(imgin,imagepathforlsd);
%imagepath='2_input.jpg';
% imagepath='img.png';
% imgout = double(imread(imagepath))/255; % This image is assumed to be the output of local warping
% Generate displacement field. 
[rows, cols, colors] = size(imgin);

% ux = zeros(rows,cols);  % the last dim is coordinate (y,x)
% uy = zeros(rows,cols);
[X,Y] = meshgrid(1:cols,1:rows);
%ux(X<200) = ux(X<200) +10;
%u = cat(3, uy, ux);


% load('xDispMap.mat');
% load('yDispMap.mat');


% xDispMap = double(xDispMap); % Chen
% yDispMap = double(yDispMap); % Chen


yDispMap = dispMap(:,:,2);
xDispMap = dispMap(:,:,1);
u = cat(3, yDispMap, xDispMap);

% load('dispMap2.mat');
% u = double(dispMap);

% warp back
% imgin = zeros(size(imgout));
% for i = 1:rows
%     for j = 1:cols
%         imgin(i+u(i,j,1), j+u(i,j,2), :)=imgout(i,j,:);
%     end
% end
% 
% imshow(imgin);

%%  Define grid mesh V
xgridN=30;
ygridN=20;
[xgrid,ygrid] = meshgrid(1:(cols-1)/(xgridN-1):cols,1:(rows-1)/(ygridN-1):rows);
xgrid = round(xgrid);
ygrid = round(ygrid);
%xgrid(10,10) = xgrid(10,10) -10;
%ygrid(10,10) = ygrid(10,10) -10;

Vl = cat(3, ygrid, xgrid);
%% Draw grid on image
gridmask = drawGridmask(Vl,  rows,cols);
imageGrided = drawGrid( gridmask, imgin);
figure(1);imshow(imageGrided);
%% Wrap grid according to displacement field
yN=size(Vl,1);
xN=size(Vl,2);
Vwarp = zeros(size(Vl));
for x =1:xN
    for y =1:yN
        Vwarp(y,x,:) = Vl(y,x,:) + u(Vl(y,x,1),Vl(y,x,2),:);
        %Vwarp(y,x,:) = u(Vl(y,x,1),Vl(y,x,2),:); % Chen
    end
end

%% Draw grid on image
gridmask = drawGridmask(Vwarp,  rows,cols);
imageGrided2 = drawGrid( gridmask, imgin);
figure(2); imshow(imageGrided2);
%imwrite(imageGrided2,[datapath 'local_mesh.png']);
%% Shape preservation matrix 

quadrows=ygridN-1;
quadcols=xgridN-1;
Vx = Vl;
Ses= cell(quadrows,quadcols);
for r = 1:quadrows
    for c = 1:quadcols
        Vo = zeros(4,2);
        Vo(1,:) = Vwarp(r,c,:); 
        Vo(2,:) = Vwarp(r,c+1,:);
        Vo(3,:) = Vwarp(r+1,c,:);
        Vo(4,:) = Vwarp(r+1,c+1,:);
        Aq = zeros(8,4);
        for i = 1:4
            Aq(i*2-1,1) = Vo(i,2);  % Vq should be x0 y0, but Vo is y x
            Aq(i*2-1,2) = - Vo(i,1);
            Aq(i*2,1) = Vo(i,1);
            Aq(i*2,2) = Vo(i,2);
            Aq(i*2-1:i*2,3:4) = eye(2);
        end
        Ses{r,c} = Aq*inv(Aq'*Aq)*Aq' - eye(8);
        
    end
end


%% boundary constraints matrix
% Vq= [x0 y0 x1 y1 ...
vertexesNum = xgridN*ygridN;
Dvec = zeros(vertexesNum*2,1);
B = zeros((vertexesNum)*2,1);
Dvec(1:xgridN*2:end) = 1;
B(1:xgridN*2:end) = 1;
Dvec(xgridN*2-1:xgridN*2:end) = 1;
B(xgridN*2-1:xgridN*2:end) = cols;
Dvec(2:2:xgridN*2) = 1;
B(2:2:xgridN*2) = 1;
Dvec(end -xgridN*2+2:2:end) = 1;
B(end -xgridN*2+2:2:end) = rows;
Dg = diag(Dvec);


%% line preservation
addpath('reference code/');
lines = lsd(imagepathforlsd);
fatmask = mask;
fatmask(1:end,1) = 1;
fatmask(1:end,end) = 1;
fatmask(1,1:end) = 1;
fatmask(end,1:end) = 1;
fatmask = filter2(ones(3,3),fatmask);
% im = imread(imagepath);
% imshow(im);
figure(1);
hold on;
for i = 1:size(lines, 2)
    plot(lines(1:2, i), lines(3:4, i), 'LineWidth', lines(5, i) / 2, 'Color', [1, 0, 0]);
end
hold off
figure(2);
% cut the lines by grid
lineN = size(lines,2);


% build warping from input to output
iu = - ones(rows,cols,2); 

for i = 1: rows
    for j = 1:cols
        iu(i+yDispMap(i,j),j+xDispMap(i,j),:) = [i j]; % Mine
        %iu(yDispMap(i,j),xDispMap(i,j),:) = [i j]; %Chen
    end
end

lineSeg = cell(quadrows, quadcols);
for i = 1:lineN
    % Warp the line points according to displace ment field.
    % then check 
    aLine=lines(1:4,i); % x1 x2 y1 y2 w
    aLine = round(aLine);
    aLine = reshape(aLine, 2, 2); % x1 y1 ; x2 y2
    aLine = fliplr(aLine); % y1 x1 ; y2 x2
    try % aline might contain negative coordinate
        if fatmask(aLine(1,1),aLine(1,2)) >0 || fatmask(aLine(2,1),aLine(2,2)) >0
            continue;
        end
    catch
        continue;
    end
    line1out = iu(aLine(1,1),aLine(1,2),:);
    line2out = iu(aLine(2,1),aLine(2,2),:);
    
    gridwidth = (cols-1)/(xgridN-1);
    gridheight = (rows-1)/(ygridN-1);
    
    line1gridID = ceil([ line1out(1)/gridheight line1out(2)/gridwidth]);
    line2gridID = ceil([ line2out(1)/gridheight line2out(2)/gridwidth]);
   
    
    hold on;
    %plot(lines(1:2,i), lines(3:4,i), 'LineWidth', lines(5,i) / 2, 'Color', [1, 0, 0]);
    hold off;
    % use getLinTrans to make sure ? 
    
    % cut lines. Follow the line from start point to end point and check
    % the boundaries it passes
    currentGrid = line1gridID;
    pstart = aLine(1,:);
    pend = aLine(2,:);
    gridstep = zeros(2,4);
    gridstep(:,1) = [0 -1];
    gridstep(:,2) = [-1 0];
    gridstep(:,3) = [1 0];
    gridstep(:,4) = [0 1];
    lineSegtmp =lineSeg;
    
    try % A line might be outside the grid, causing error. Such line should be discarded.
        lineplot1 = 0;
        lineplot2 = 0;
        while(true)
            
            if currentGrid == line2gridID
                pend = aLine(2,:);
            else

                Vq=getvertices(currentGrid, Vwarp);

                quadlines = zeros(2,2,4);
                quadlines(:,:,1) = [Vq(1,:); Vq(3,:)];
                quadlines(:,:,2) =  Vq(1:2,:);
                quadlines(:,:,3) = Vq(3:4,:);
                quadlines(:,:,4) = [Vq(2,:); Vq(4,:)];
                findInterSect=false;
                for l = 1:4
                    [isInter, p]= checkInterSection(quadlines(:,:,l), [pstart ; aLine(2,:)]);
                    if isInter  && norm(pstart-p) >0
                        %norm(pstart-p)
                        pend = p;
                        hold on;
                        lineplot1 = plot(quadlines(:,2,l),quadlines(:,1,l), 'LineWidth', 1, 'Color', [0, 0, 1]);
                        hold off;
                        nextGrid = currentGrid +gridstep(:,l)';
                        findInterSect = true;
                    end
                end
                if (findInterSect==false)
                    break;
                end
            end
            
            %% add lines to quads
            hold on;
            lineplot2 = plot([pstart(2) pend(2)], [pstart(1) pend(1)], 'LineWidth', lines(5,i) / 2, 'Color', [1, 0, 0]);
            % some lines are missing
            hold off;
            lineSegtmp{currentGrid(1),currentGrid(2)}= [lineSegtmp{currentGrid(1),currentGrid(2)} ; [pstart pend 0]]; % 0 is the rotation angle
            
            if (currentGrid == line2gridID)
                break;
            end
            currentGrid = nextGrid;
            pstart = pend;
            if norm(pstart-aLine(2,:)) <1;
                break;
            end
        end
        
    catch 
        if lineplot1
            delete(lineplot1);
        end
        if lineplot2
            delete(lineplot2);
        end
        continue;
    end
    lineSeg =lineSegtmp;
    
end

%% Build shape matrix
% vertexes to quads matrix
quadrows=ygridN-1;
quadcols=xgridN-1;
Q = zeros(8*quadrows*quadcols,2*ygridN*xgridN);
for i = 1:quadrows
    for j = 1:quadcols
        quadID = ((i-1)*quadcols+j-1)*8;
        topleftvertexID = ((i-1)*xgridN+j-1)*2;
        Q(quadID+1:quadID+2, topleftvertexID+1:topleftvertexID+2) = eye(2); % x0 y0 topleft vertex
        Q(quadID+3:quadID+4, topleftvertexID+3:topleftvertexID+4) = eye(2); % x1 y1 topright vertex
        Q(quadID+5:quadID+6, topleftvertexID+xgridN*2+1:topleftvertexID+xgridN*2+2) = eye(2); % x2 y2 bottom left vertex
        Q(quadID+7:quadID+8, topleftvertexID+xgridN*2+3:topleftvertexID+xgridN*2+4) = eye(2); % x3 y3 bottom right vertex
    end
end
S = []; % shape
for i = 1:quadrows
    Si = [];
    for j = 1:quadcols
        Si = blkdiag(Si,Ses{i,j});
    end
    S = blkdiag(S,Si);
end

%% calculate theta of original line
% group line segments
linegroup = cell(quadrows, quadcols);
qstep = pi/49;
for i = 1:quadrows
    for j = 1:quadcols
        lineNum = size(lineSeg{i,j},1);
        for l = 1:lineNum
            pstart = lineSeg{i,j}(l,2:-1:1)'; % x y
            pend = lineSeg{i,j}(l,4:-1:3)';
            theta = atan((pstart(2)-pend(2))/(pstart(1)-pend(1)));
            groupid = round((theta+pi/2)/qstep)+1;
            linegroup{i,j} = [linegroup{i,j}; [groupid theta]];
        end
    end
end


%% optimization loop
itNum = 10;
for it = 1: itNum
    
    %% iterate over all lines in quads to create the C matrixes  (recompute in each iteration)
    % quantize orienation range
    Cmatrixes = cell(quadrows, quadcols);
    LinTrans = cell(quadrows, quadcols);
    for i = 1:quadrows
        for j = 1:quadcols
            lineNum = size(lineSeg{i,j},1);
            for l = 1:lineNum
                Vq=getvertices([i j], Vwarp); %This Vq is y x
                Vq = fliplr(Vq);
                Vq = reshape(Vq',[],1); % Vq should be x y
                pstart = lineSeg{i,j}(l,2:-1:1)';
                pend = lineSeg{i,j}(l,4:-1:3)';
                T1 = getLinTrans( Vq,pstart); % lineSeg is y x
                T2 = getLinTrans( Vq,pend );
                LinTrans{i,j} = [LinTrans{i,j}; [T1 T2]];
                ehat = pstart - pend;
                theta = lineSeg{i,j}(l,5);
                R = [cos(theta) -sin(theta) ; sin(theta) cos(theta)];
                C= R*ehat*inv(ehat'*ehat)*ehat'*R' -eye(2);
                CT = C*(T1-T2);
                Cmatrixes{i,j} = [Cmatrixes{i,j} ; CT];
            end
        end
    end
    
    
    %% update line matrix includes H (quads to lines)
    L = [];
    Nl = 0;
    for i = 1:quadrows
        Li =[];
        for j = 1:quadcols
            lineNum = size(lineSeg{i,j},1);
            Nl = Nl +lineNum;
            if lineNum == 0
                Li = [Li zeros(size(Li,1),8)];
            else
                Li = blkdiag(Li,Cmatrixes{i,j});
            end
        end
        L = blkdiag(L,Li);
    end
    
    %% combine matrixes
    Nq = quadrows*quadcols;
    lambL = 100;
    lambB = 1e8;
    K = [ 1/Nq*S*Q;  lambL/Nl*L*Q; lambB*Dg];
    BA = [zeros(size(K,1) - size(B,1),1); lambB*B];
    
    %% Update V solving linear system
    
    A = K'*K;
    b = K'*BA;
    x = A\b;
    Vopt = zeros(ygridN,  xgridN,2);
    %% draw
    for i = 1:ygridN
        for j = 1:xgridN
            xid = ((i-1)*xgridN + j-1)*2;
            Vopt(i,j,:) = x(xid+2:-1:xid+1);
        end
    end
    
    %gridmask = drawGridmask(Vopt,  rows,cols);  %% FIXME
    %imageGrided2 = drawGrid( gridmask, imgin);
    %figure(3); imshow(imageGrided2);

    %%  calculate the angle of new line segs
    thetagroup = zeros(50,1);
    thetagroupNum = zeros(50,1);
    for i = 1:quadrows
        for j = 1:quadcols
            lineNum = size(lineSeg{i,j},1);
            Vq=getvertices([i j], Vopt); %This Vq is y x
            Vq = fliplr(Vq);
            Vq = reshape(Vq',[],1); % Vq should be x y
            for l = 1:lineNum
                T1 = LinTrans{i,j}((l-1)*2+1:(l-1)*2+2,1:8);
                T2 = LinTrans{i,j}((l-1)*2+1:(l-1)*2+2,9:16);
                pstartnew = T1*Vq;
                pendnew = T2*Vq;
                theta = atan((pstartnew(2)-pendnew(2))/(pstartnew(1)-pendnew(1)));
                deltatheta = theta-linegroup{i,j}(l,2);
                if isnan(linegroup{i,j}(l,2)) || isnan(deltatheta)
                    continue;
                end
                if deltatheta > pi/2
                    deltatheta = deltatheta - pi;
                end
                if deltatheta < -pi/2
                    deltatheta = pi + deltatheta;
                end
                
                thetagroup( linegroup{i,j}(l,1)) = thetagroup( linegroup{i,j}(l,1))+deltatheta;
                thetagroupNum( linegroup{i,j}(l,1)) = thetagroupNum( linegroup{i,j}(l,1)) +1;
                % plot new lines
                %             hold on;
                %                 lineplot2 = plot([pstartnew(1) pendnew(1)], [pstartnew(2) pendnew(2)], 'Color', [1, 0, 0]);
                %             hold off;
                
            end
        end
    end
    %% calculate mean theta of each bin
    thetagroup = thetagroup./thetagroupNum;
    for i = 1:quadrows
        for j = 1:quadcols
            lineNum = size(lineSeg{i,j},1);
            for l = 1:lineNum
                lineSeg{i,j}(l,5) = thetagroup(linegroup{i,j}(l,1));
            end
        end
    end

end  % optimization iteration end


% %% perform warping based on Vwarp and Vopt 
% outputimg = zeros(rows,cols,3);
% outputimgcount = zeros(rows,cols,3);
% currentquadID = [1 1];
% nextquadID = [1 1];
% % for i = 1:rows
% %     for j = 1:cols
% %         % locate quad
% %         Vq=getvertices(currentquadID, Vopt);
% %         Vq = fliplr(Vq);
% %         Vq = reshape(Vq',[],1); % Vq should be x y
% %         Tp = getLinTrans( Vq,[j i]');
% %         t1n = Tp(1,5)+Tp(1,7);
% %         t2n = Tp(1,3)+Tp(1,7);
% %         if t1n< -1e-8
% %             nextquadID = currentquadID + [-1 0];
% %         end
% %         if t1n > 1 +1e-8
% %             nextquadID = currentquadID + [1 0];
% %         end
% %         if t2n< -1e-8 
% %             nextquadID = currentquadID + [0 -1];
% %         end
% %         if t2n> 1 +1e-8
% %             nextquadID = currentquadID + [0 1];
% %         end   
% %         if ~ isequal(currentquadID,nextquadID)
% %             j = j -1;
% %             currentquadID = nextquadID;
% %             figure(4); imshow(outputimg);
% %             continue;
% %         end
% %         
% %         Vo=getvertices(currentquadID, Vwarp);
% %         Vo = fliplr(Vo);
% %         Vo = reshape(Vo',[],1); % Vq should be x y
% %         pref = round(Tp*Vo); 
% %         
% %         outputimg(i,j,:)=imgout(pref(2),pref(1),:);
% %         if j == cols
% %             currentquadID(2) =1; 
% %         end
% %     end
% % end
% 
% for i = 1:quadrows
%     for j = 1:quadcols
%         Vq=getvertices([i j], Vopt);
%         Vq = fliplr(Vq);
%         Vq = reshape(Vq',[],1); % Vq should be x y
%         Vo=getvertices([i j], Vwarp);
%         Vo = fliplr(Vo);
%         Vo = reshape(Vo',[],1); % Vq should be x y
%         xlength = max(Vq(1:2:8)) -min(Vq(1:2:8));
%         ylength = max(Vq(2:2:8)) -min(Vq(2:2:8));
%         t2nstep = 1/(xlength*4);
%         t1nstep = 1/(ylength*4);
%         for t1n = 0:t1nstep:1
%             for t2n = 0:t2nstep:1
%                 v1w=1-t1n-t2n+t1n*t2n;
%                 v2w=t2n-t1n*t2n;
%                 v3w=t1n-t1n*t2n;
%                 v4w=t1n*t2n;
%                 T = [ v1w 0 v2w 0 v3w 0 v4w 0; ...
%                         0 v1w 0 v2w 0 v3w 0 v4w ];
%                 pout = round(T*Vq);
%                 pref = round(T*Vo); 
%                 try  % Vq may contain negetive index
%                 outputimg(pout(2),pout(1),:) = outputimg(pout(2),pout(1),:)+ imgin(pref(2),pref(1),:);
%                 outputimgcount (pout(2),pout(1),:) = outputimgcount (pout(2),pout(1),:) +1;
%                 catch
%                     
%                 end
%             end
%         end
%         figure(4); imshow(outputimg./outputimgcount);
%     end
% end
% outputimg= outputimg./outputimgcount;
% imwrite(outputimg,[datapath 'global.png']);
% 
% gridmask = drawGridmask(Vopt,  rows,cols);
% imageGrided2 = drawGrid( gridmask, outputimg);
% figure(5); imshow(imageGrided2);
% imwrite(imageGrided2,[datapath 'global_mesh.png']);

Vlocal = Vwarp;
Vglobal = Vopt;
end

