function str_exit = FEMcheeger2(nphi,ff,alpha,ppn,form,start_point)

% --- June 15th 2017 --------------
% Algorithm based on the article
% "Phase Field Approach to Optimal Packing Problems and 
% Related Cheeger Clusters"
% by Beniamin Bogosel, Dorin Bucur and Ilaria Fragala
% 
% Syntax: 
% FEMcheeger2(nphi,ff,alpha,ppn,form,start_point)
%     nphi - number of cells in the Cheeger cluster. 
%            chosse nphi=1 if you want to compute a Cheeger cell
%     ff   - epsilon parameter in Gamma-convergence. Generally ff=1 works well
%     alpha - parameter in the definition of alpha Cheeger.
%             alpha=1 corresponds to the classical Cheeger set.
%     ppn   - parameter for the p-norm. Choose 1 for the sum, 100 for the maximum
%     form  - 'equi', 'square' load predefined meshes
%           - otherwise, a mesh is generated for the polygon 
%             defined by pts below
%     start_point - optional. Give a structure containing a starting point
%                 - useful to continue or refine an existing optimization
%
% Examples:
% Cheeger set for equilateral triangle
%    FEMcheeger2(1,1,1,1,'equi')   
%
% 0.5001-Cheeger set for the square - should give something close to the incircle
%    FEMcheeger2(1,1,0.5001,1,'square')
% 
% 10-Cheeger cluster for equilateral triangle
%    FEMcheeger2(10,1,1,1,'equi')   
%
% 10-circle packing for equilateral triangle
%    FEMcheeger2(10,1,0.5001,1,'equi')   
%
%  depends on  - LBFGS or minConf for the optimization 
%              - mesh2D if you want to mesh general domains
%              - tricontour if you want to plot contour in the end

tic
clf

method = 2;  % choose optimization method
             % 1 - LBFGS 
             % 2 - minConf (faster)
fact = 100;  % penalization parameter for the phase separation


switch form
case 'equi'
 load mesh_equi;
 mesh_struc = mesh_str;
 pts = [1 0; cos(2*pi/3) sin(2*pi/3); cos(4*pi/3) sin(4*pi/3)];
 ai = 3*sqrt(3)/4;
 pts  = sqrt(pi/ai)*pts;
case 'square'
 load mesh_square;
 mesh_struc = mesh_str;
 pts = [-1 -1;
       1 -1;
       1 1;
       -1 1];
otherwise
 nphi = 1;
 %%%%%%%% modify here if you want some particular domain in mind %%%%%%%%%%
 pts = 2*[-0.6 -0.3;
               1 0;
               1.7 1;
               0.5 2;
              -1 0;];
 %pts = [0 0; 1 0; 1 1; 0 1]; %-0.3 0.5];
 %pts = 3*[1 0; cos(2*pi/3) sin(2*pi/3); cos(4*pi/3) sin(4*pi/3)]; 
 [p,t] = mesh2D_poly(pts);
 ms.points = p;
 ms.t      = t;
 ms.dx     = dx_calculator(p',t);
 %pause
 ms.Ibord  = Ibord_calc(p,t);
 [K,M]             = dir_assemKM(p,t);
 ms.K      = K;
 ms.M      = M;
 ms.npt    = size(p,1);
 mesh_struc{1}     = ms;
end

nref = size(mesh_struc,2);


for i = 1:nref
if i>1 
 ff=1;
end

if i>0
 pn = ppn;
else
 pn = 1;
end


astruc = mesh_struc{i}

points  = astruc.points;
npt     = astruc.npt
Prob.K  = astruc.K;
Prob.M  = astruc.M;
Prob.ep = ff*astruc.dx;
dx      = astruc.dx;
Prob.points = points;
Prob.nphi   = nphi;
Prob.npt    = npt;
Prob.alpha  = alpha;
Prob.fact   = fact;
Prob.pn     = pn;
Prob.Ibord  = astruc.Ibord;
Ibord       = astruc.Ibord;

Prob

% admissible starting point

if i==1
 if(exist('start_point','var')==0)
   M0                                  = rand(npt,nphi);
 else
   M0                                  = start_point.U;
   M0 = reshape(M0,npt,nphi);
 end
else
   ref = astruc.ref;
   M0 = U(ref,:);
end
 if i>1
  %M0 = sparse(M0);
 end 
  M0(Ibord,:) = 0;
  vec0 = M0(:);
  Prob.optParam.MaxIter                = 5000;
  opts                                 = lbfgs_options('iprint', -1, 'maxits', Prob.optParam.MaxIter , 'cb', @spartitions_callback,'pgtol',1e-10,'m',5);
  
  lbd                                   = zeros(size(vec0));
  ubd                                  = ones(size(vec0));
  nbd                                  = 2*ones(size(vec0));
  
  opts.Prob                            = Prob;
  opts.Prob.factaway = 0;
  
  disp('Test derivative:\n');

  u           = M0(:);
  direction   = rand(size(u));
  ep          = 1e-6;

  [val,gradt] = FEMch_cost(u,Prob);

  val1        = FEMch_cost(u+ep*direction,Prob);
  val2        = FEMch_cost(u-ep*direction,Prob);

  %fprintf('Analytical gradient: %6.6f\n',sum(gradt.*direction));
  %fprintf('Finite differences:  %6.6f\n',(val1-val2)/(2*ep));
  %pause
decrease = 1;




switch method
case 1
   [ U,fx,exitflag,userdata]  = LBFGS(@FEMch_cost,vec0,lbd,ubd,nbd,opts);
U = reshape(U,[],nphi);
case 2
   options.method = 'lbfgs';
   options.optTol = 1e-15;
   %options.maxIter = 1000;
   [U,fx] = minConf_TMP(@(x) FEMch_cost(x,Prob),vec0,lbd,ubd,options);
   U = reshape(U,[],nphi);
end

  
end

U = U(:);

str_exit.U         = U;
str_exit.points    = points;
str_exit.val       = fx;
str_exit.t         = astruc.t;
str_exit.dx        = dx;
str_exit.pts       = pts;
p = astruc.points';


y = reshape(U,[],nphi);
if nphi>1
[I,J] = max(y');
else
I   = U;
end
size(astruc.points)
figure(1)
clf
patch('Faces',astruc.t,'Vertices',astruc.points,'FaceVertexCData',I(:),'FaceColor','interp','EdgeColor','none');

view(2)
axis equal
axis off
toc
   size(I)
   lev = 0.2;
   tricontour(astruc.points,astruc.t,I(:),[lev lev]);
if nphi ==1
  try
   figure(2)
   clf
   pts = [pts;
          pts(1,:)];
   plot(pts(:,1),pts(:,2),'LineWidth',2);
   lev = dx;
   tricontour(astruc.points,astruc.t,I,[lev lev],3)
   ipts = cheeger_poly(pts);
   plot(ipts(:,1),ipts(:,2),'--','LineWidth',4);
   hold off
  end
end


