function ptch = patchSphere(sfit,N)
% PATCHSPHERE creates a triangle-based patch of a sphere with 
% equidistributed points on the surface of the sphere.
%
%   ptch = patchSphere(sfit) uses a default number of points approximately
%   equal to 100 (may vary due to rounding).
%
%       sfit - a sphereModel object or one of the following representations
%              of a sphere 
%               (1) a 1x1 array containing the radius of the sphere only; 
%               (2) a 1x4 array containing sphere parameters [a,b,c,r] such
%                   that (x-a)^2 + (y-b)^2 + (z-c)^2 = r^2; OR
%               (3) a structured array containing the fields "Center" and
%                   "Radius".
%
%   ptch = patchSphere(sfit,N) uses a default number of points 
%   approximately equal to N (may vary due to rounding.
%
%   References:
%       [1] Bentz, B. "mySphere," 2016.
%       https://www.mathworks.com/matlabcentral/fileexchange/57877-mysphere-n
%
%       [2] Deserno, M. "How to generate equidistributed points on the
%       surface of a sphere", 2004.
%       https://www.cmu.edu/biolphys/deserno/pdf/sphere_equi.pdf
%
%   M. Kutzer, 15Jun2020, USNA

%% Check inputs and default value(s)
narginchk(1,2);

% Standardize input
sfit = sphereParam2Model(sfit);

if nargin < 2
    N = 100;
end

%% Parse sphere radius
r = sfit.Radius;

%% The following code is the code of B. Bentz with one tiny changes
r_unit = 1;

Area = 4*pi*r_unit^2/N;
Distance = sqrt(Area);
M_theta = round(pi/Distance);
d_theta = pi/M_theta;
d_phi = Area/d_theta;

N_new = 0;
for m = 0:M_theta-1
    
    Theta = pi*(m+0.5)/M_theta;
    M_phi = round(2*pi*sin(Theta)/d_phi); % not exact
    
    for n = 0:M_phi-1        
        Phi = 2*pi*n/M_phi;    
        
        N_new = N_new + 1;
        
        X(N_new) = r*sin(Theta)*cos(Phi);
        Y(N_new) = r*sin(Theta)*sin(Phi);
        Z(N_new) = r*cos(Theta);
        
    end
end

% figure;
% plot3(X.',Y.',Z.','.m');

%% Create triangle-based mesh
DT = delaunayTriangulation(X.',Y.',Z.');
[K,~] = convexHull(DT);
ptch.Faces = K;
ptch.Vertices = DT.Points;
%patch(ptch,'FaceColor','b','EdgeColor','k','FaceAlpha',0.5);

%% Center sphere
X0 = ptch.Vertices;
X = X0 + reshape(sfit.Center,1,[]);
ptch.Vertices = X;
