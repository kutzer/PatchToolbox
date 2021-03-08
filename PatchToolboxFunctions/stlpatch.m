function [p,stltitle] = stlpatch(filename)
% stlpatch creates a structured variable formatted as for patch. 
%   stlpatch prompts the user for a binary STL file. The binary STL is
%   read and formatted as a "patch." STL title information is also
%   returned.
%
%   stlpatch(filename) reads the binary STL specified using filename and
%   creates a structured variable formatted as a "patch." STL title
%   information is also returned.
%
%   See also patch, stlread
%
%   M. Kutzer, 19Dec2014, USNA

%TODO - check for stlread and point to website:
% http://www.mathworks.com/matlabcentral/fileexchange/29906-binary-stl-file-reader

[v, f, ~, c, stltitle] = stlread(filename);

p.FaceColor = mean(c); % patch objects do not currently enable individual face colors
p.FaceAlpha = 1;
p.EdgeColor = 'none';
p.LineStyle = 'none';
p.Faces = f;
p.Vertices = v;