function [p,t] = mesh2D_poly(pts)

if nargin<1
  pts = [-1 -0.1;
         1 0;
         1 1;
         0.5 2;
         -1 1;];
end


npts = size(pts,1);

ind = (1:npts)';
ind2 = circshift(ind,-1);
edge = [ind ind2]

node = pts;
hfun = 0.1;

opts.kind = 'delfront';
opts.rho2 = +1 ;
%opts.dtri = 'conforming';
%opts.siz2 = 0.01 ;
%opts.disp = Inf;

[vert,etri, ...
    tria,tnum] = refine2(node,edge,[],opts,0.02) ;
[vert,etri, ...
    tria,tnum] = smooth2(vert,etri,tria,tnum) ;
    

size(vert)
size(etri)
size(tria)

clf
patch('faces',tria(:,1:3),'vertices',vert, ...
        'facecolor','w', ...
        'edgecolor',[.2,.2,.2]) ;
axis equal

p = vert;
t = tria(:,1:3);
