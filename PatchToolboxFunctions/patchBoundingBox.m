function p = patchBoundingBox(X_lims)
% PATCHBOUNDINGBOX visualizes the bounding box for a given set of limits
%   p = PATCHBOUNDINGBOX(X_lims) visualizes the bounding box described by
%   X_lims.
%
%   Inputs:
%       X_lims - 3x2 array containing x/y/z limits
%   Outputs:
%       p - structured array containing patch fields
%
%   M. Kutzer & L. Davis, 04Feb2021, USNA

%% Debug plot
debugON = false;
if debugON
    fig = figure;
    axs = axes('Parent',fig);
    hold(axs,'on');
    daspect(axs,[1 1 1]);
    view(axs,3);
end

%% Check input(s)
narginchk(1,1);
if size(X_lims,1) ~= 3 || size(X_lims,2) ~= 2
    error('The input must be a 3x2 array containing x/y/z limits.');
end

%% Define vertices
x = X_lims(1,:);
y = X_lims(2,:);
z = X_lims(3,:);

v = [...
    x(1),y(1),z(1);...
    x(2),y(1),z(1);...
    x(2),y(2),z(1);...
    x(1),y(2),z(1);...
    x(1),y(1),z(2);...
    x(2),y(1),z(2);...
    x(2),y(2),z(2);...
    x(1),y(2),z(2)];

if debugON
    plt = plot3(axs,v(:,1),v(:,2),v(:,3),'*k');
    for i = 1:size(v,1)
        txt(i) = text(v(i,1),v(i,2),v(i,3),sprintf('%d',i),...
            'HorizontalAlignment','Left','VerticalAlignment','Top',...
            'Parent',axs);
    end
end

%% Define Faces
f = [
    1,2,6,5;... % 1
    5,6,7,8;... % 2
    1,5,8,4;... % 3
    2,3,7,6;... % 4
    1,4,3,2;... % 5
    3,4,8,7];   % 6

if debugON
    colors = 'rgbcmy';
    for i = 1:size(f,1)
        % Plot Faces
        ptc(i) = patch('Parent',axs,'Faces',f(i,:),'Vertices',v,...
            'FaceColor',colors(i),'EdgeColor','none','FaceAlpha',0.5);
        m = mean( v(f(i,:),:), 1 );
        txtc(i) = text(axs,m(1),m(2),m(3),sprintf('f_%d',i));
        % Plot normals
        v1 = v(f(i,2),:) - v(f(i,1),:);
        v2 = v(f(i,4),:) - v(f(i,1),:);
        n = cross(v1,v2);
        n_hat = n./norm(n);
        qvr(i) = quiver3(m(1),m(2),m(3),n_hat(1),n_hat(2),n_hat(3),...
            'Parent',axs,'Color',colors(i),'LineWidth',2);
    end
end
    
%% Package patch struct
p.Vertices = v;
p.Faces = f;
    