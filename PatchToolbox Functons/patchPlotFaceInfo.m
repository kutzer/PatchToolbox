function txt = patchPlotFaceInfo(ptch)
% PATCHPLOTFACEINFO overlays the individual face number and vertex
% order for each face of a provided patch object.
%   txt = PATCHPLOTFACEINFO(ptch) overlays the face number and vertex order
%   for a patch object ptch. An array of text handles is returned. 
%
%   M. Kutzer, 07May2019, USNA

% Get parent of patch object 
mom = get(ptch,'Parent');

% Get vertices and faces of patch object
v = ptch.Vertices;
f = ptch.Faces;

% Create temporary axes
figTMP = figure('Visible','off');
axsTMP = axes('Parent',figTMP);
hold(axsTMP,'on');

for i = 1:size(f,1)
    for j = 1:size(f,2)
        f_now = f(i,j);
        v_now = v(f_now,:);
        
        str = sprintf('%d_{%d}',i,j);
        
        txt(i,j) = text(v_now(1),v_now(2),v_now(3),str,'Parent',axsTMP);
    end
end

set(txt,'Parent',mom);
delete(figTMP);