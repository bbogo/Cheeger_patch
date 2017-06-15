function ipts = cheeger_poly(pts,arg2)

% use Lachand-Robert - Kawohl formulation to find the Cheeger set
% !! find t such that the area of the t-inner offset = pi*t^2
%
% relies on the Clipper Library for computing efficiently polygon offsets
% 
% Syntax
%    cheeger_poly(pts,arg2)
%        pts - 2xN matrix containing coordinates of points in 2D
%        arg2- additional argument. If used, a plot of the polygon and
%              of the Cheeger set appears

if nargin>1
 plotting = 1;
else
 plotting = 0;
end

if nargin<1
 pts = [-1 0.6;
         1 0;
         1 1;
         0.5 2;
        -1 1;];
end

%dis = linspace(0,1,100);
%ys  = arrayfun(@(x) cheeger_diff(pts,x),dis);

%if plotting==1
 %figure(2)
 %plot(dis,ys);
%end

% search for the zero of cheeger_diff by bissection

x = fzero(@(t) cheeger_diff(pts,t),0);
%if plotting==1
% hold on
% plot(x,0,'or');
%end
qpts = [pts; pts(1,:)];
[X,Y] = polyout(pts(:,1),pts(:,2),-x);
ipts = [X{:} Y{:}];
ipts = [ipts; ipts(1,:)];


if plotting==1
  clf
  plot(qpts(:,1),qpts(:,2),'k','LineWidth',2);
  hold on
  plot(ipts(:,1),ipts(:,2),'r','LineWidth',2);
end
  [X,Y] = polyout(ipts(:,1),ipts(:,2),x,'r',0.0001);
  ipts = [X{:} Y{:}];
  ipts = [ipts; ipts(1,:)];
if plotting==1
  plot(ipts(:,1),ipts(:,2),'LineWidth',2);
  axis equal
end

fprintf('Cheeger constant = %.10f\n',1/x);
fprintf('Inverse          = %.10f\n',x);

function res = cheeger_diff(pts,t)

[X,Y] = polyout(pts(:,1),pts(:,2),-t);
if prod(size(X))>0
ar    = polyarea(X{:},Y{:});
res   = ar-pi*t^2;
else
res   = -pi*t^2;
end
