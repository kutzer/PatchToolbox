function varargout = patchBoundingBox(varargin)
% PATCHBOUNDINGBOX defines face and vertex information for a bounding box 
% given set of limits OR defines a bounding box given one or more patch
% objects.
%   p = PATCHBOUNDINGBOX(X_lims) defines face and vertex information for a 
%   bounding box given a set of limits described by X_lims.
%
%   Inputs:
%       X_lims - 3x2 array containing x/y/z limits
%   Outputs:
%       p - structured array containing patch fields
%
%   [X_lims, X_lim_i] = PATCHBOUNDINGBOX(p) calculates the overall bounding
%   box for one or more patch objects (X_lims), and calculates the
%   individual bounding box for each element of p (X_lim_i).
%
%   Inputs:
%       p - n-element structured array containing fields "Faces" and
%           "Vertices" or a valid array of patch objects.
%
%   Outputs:
%       X_lims  - 3x2 array containing x/y/z limits for all patch objects
%       X_lim_i - n-element cell array
%           X_lim_i{i} - 3x2 array containing x/y/z limits for the ith
%                        patch object
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
switch class(varargin{1})
    case 'struct'
        p = varargin{1};
        X_lims = [];
        fields = {'Faces','Vertices'};
        bin = isfield(p,fields);
        if nnz(bin) ~= numel(fields)
            error('Input must be a struct with fields "Faces" and "Vertices" or a valid patch object');
        end
    case 'matlab.graphics.primitive.Patch'
        p = varargin{1};
        X_lims = [];
    otherwise
        X_lims = varargin{1};
end
 
%% Define X_lims (if patch is provided)
if isempty(X_lims)
    for i = 1:numel(p)
        v = p(i).Vertices;
        for j = 1:size(v,2)
            X_lim_i{i}(j,:) = [min(v(:,j)),max(v(:,j))];
        end
        
        X_lims = [X_lims, X_lim_i{i}];
    end
    X_lims = [min(X_lims,[],2),max(X_lims,[],2)];
    
    % Package outputs
    varargout{1} = X_lims;
    varargout{2} = X_lim_i;
    return
end

%% Run remainder of code
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

varargout{1} = p;
    