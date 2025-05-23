%% SCRIPT_Example_RotateTimberWolf
clear all
close all
clc

%% Load 3D visualization
% Open figure
fig = open('Timberwolf.fig');
% Make figure visible
set(fig,'Visible','on');
% Recover axes
axs = get(fig,'Children');
% Recover transform
h_o2w = findobj('Parent',axs,'Type','hgtransform');
% Recover patch
ptc = findobj('Parent',axs,'Type','patch');

%% Adjust view & axes limits
view(axs,[60,30]);
xlim(axs,[-2,2]);
ylim(axs,[-2,2]);
zlim(axs,[-2,2]);

%% Adjust axes size in figure
set(axs,'Units','Normalized');
set(axs,'Position',[0,0,1,1]);

%% Adjust parent/child relationship(s)
set(ptc,'Parent',h_o2w);

%% Make a video of the robot rotating
fps = 30;
vidDuration = 3; % seconds
thetas = linspace(0,2*pi, round(fps * vidDuration));
%thetas(end) = [];   % Remove the last angle for a continuous rotating gif

% Setup the video object
vidObj = VideoWriter('Timberwolf_Rotate.mp4','MPEG-4');
open(vidObj);

% Remove ugly gray
set(fig,'Color',[1,1,1]);
set(axs,'Visible','off');

% Make figure big
set(fig,'Units','Normalized','Position',[0,0,0.80,0.80]);
centerfig(fig);

for theta = thetas
    set(h_o2w,'Matrix',Rx(theta));
    drawnow;
    
    % Grab a frame for the video
    frame = getframe(fig);
    % Write the frame to the video
    writeVideo(vidObj,frame);
end
for theta = thetas
    set(h_o2w,'Matrix',Ry(theta));
    drawnow;
    
    % Grab a frame for the video
    frame = getframe(fig);
    % Write the frame to the video
    writeVideo(vidObj,frame);
end
for theta = thetas
    set(h_o2w,'Matrix',Rz(theta));
    drawnow;
    
    % Grab a frame for the video
    frame = getframe(fig);
    % Write the frame to the video
    writeVideo(vidObj,frame);
end
for theta = thetas
    set(h_o2w,'Matrix',Rx(theta)*Ry(theta)*Rz(theta));
    drawnow;
    
    % Grab a frame for the video
    frame = getframe(fig);
    % Write the frame to the video
    writeVideo(vidObj,frame);
end
% Close the video object
close(vidObj);
