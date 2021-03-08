function [p,p_bb] = centerAndPatchSTL(stlFile)
% CENTERANDPATCHSTL centers an stl in x/y and creates a patch struct
%   [p,p_bb] = CENTERANDPATCHSTL(stlFile)
%
%   Inputs
%       stlFile - string with file information for stl
%   Outputs
%       p - patch struct for stl
%       p_bb - patch struct for bounding box
%
%   M. Kutzer & L. Davis, 22Feb2021, USNA

debugON = false;

%% Create patch from STL
p = stlpatch(stlFile);

%% Find the bounding box of the STL
v = p.Vertices;
for i = 1:size(v,2)
    X_lims(i,:) = [min( v(:,i) ), max( v(:,i) )];
end

p_bb = patchBoundingBox(X_lims);

%% Define a transformation to center the patch at [0,0] in [x,y]
H = Tx( X_lims(1,1) + diff(X_lims(1,:))/2 ) *...
    Ty( X_lims(2,1) + diff(X_lims(2,:))/2 ) *...
    Tz( X_lims(3,1));

v = v.';    % Make the vertices a 3xN
v(4,:) = 1; % Make the vertices homogeneous

v = H^(-1) * v;

%% Package updated patch vertices
p_tmp = p;  % Save original
p.Vertices = v(1:3,:).'; % Updated file

%% Debug code
if debugON
    fig = figure;
    axs = axes('Parent',fig);
    hold(axs,'on');
    view(axs,3);
    daspect(axs,[1 1 1]);
    addSingleLight(axs);
    
    ptc = patch(p,'FaceColor','b','EdgeColor','none','FaceAlpha',0.5);
    X_lims = [xlim(axs); ylim(axs); zlim(axs)];
    d = diff(X_lims.').';
    hg = triad('Parent',axs,'Scale',0.8*min(d),'LineWidth',2);
    
    ptc_tmp = patch(p_tmp,'FaceColor','g','EdgeColor','none','FaceAlpha',0.5);
    
    ptc_bb = patch(p_bb,'FaceColor','k','EdgeColor','k','FaceAlpha',0.1);
    hg = triad('Parent',axs,'Scale',0.8*min(d),'LineWidth',2,'Matrix',H);
end