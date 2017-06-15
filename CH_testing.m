function CH_testing

load('test_caseCH.mat','square');
whos square

% square

figure(1)
% first argument of cheeger_poly: matrix containing points: Nx2
% if more than one argument then plotting of the solution occurs

% as output arguments cheeger_poly returns the vertices of an approximation
% of the Cheeger set, to be able to plot things afterwards

P = square;
cheeger_poly(P,1);

% poly1
figure(2)
load('test_caseCH.mat','poly1');
P = poly1;
cheeger_poly(P,1);

% poly2
figure(3)
load('test_caseCH.mat','poly2');
P = poly2;
cheeger_poly(P,1);

% Reuleaux triangle
% for non-polygonal convex sets a sufficiently fine polygonal approximation 
% can be considered
% in this case the Reuleaux triangle is approximated by a polygon with 402 vertices

figure(4)
load('test_caseCH.mat','pts');
P = pts;
cheeger_poly(P,1);
