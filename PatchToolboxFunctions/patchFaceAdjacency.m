function [adj,g] = patchFaceAdjacency(ptch)
% PATCHFACEADJACENCY returns the adjacency matrix that identifies
% connections between faces with shared edges. 
%   adj = PATCHFACEADJACENCY(ptch) defines the patch as either a structured
%   array containing fields "Faces" and "Vertices" or a patch object.
%
%   NOTE: This function assumes the patch is defined using a triangular
%   mesh
%
%   M. Kutzer, 01Jun2020, USNA