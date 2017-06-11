function str_exit = FEMcheeger2(nphi,ff,alpha,ppn,form,start_point)
tic
clf
fact = 100;

switch form
case 'equi'
 load mesh_equi;
 mesh_struc = mesh_str;
case 'square'
 load mesh_square;
 mesh_struc = mesh_str;
end

nref = size(mesh_struc,2);

for i = 1:nref
if i>1 
 ff=1;
end

if i>2
 pn = ppn;
else
 pn = 1;
end

astruc = mesh_struc{i}
points  = astruc.points;
npt     = astruc.npt
Prob.K  = astruc.K;
Prob.M  = astruc.M;
Prob.ep = ff*astruc.dx
Prob.points = points;
Prob.nphi   = nphi;
Prob.npt    = npt;
Prob.alpha  = alpha;
Prob.fact   = fact;
Prob.pn     = pn;
Prob.Ibord  = astruc.Ibord;
Ibord       = astruc.Ibord;



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

  fprintf('Analytical gradient: %6.6f\n',sum(gradt.*direction));
  fprintf('Finite differences:  %6.6f\n',(val1-val2)/(2*ep));
  %pause
decrease = 1;

[ U,fx,exitflag,userdata]             = LBFGS(@FEMch_cost,vec0,lbd,ubd,nbd,opts);
U = reshape(U,[],nphi);
  
end

U = U(:);

str_exit.U         = U;
str_exit.points    = points;
str_exit.val       = fx;
str_exit.t         = astruc.t;
p = astruc.points';


y = reshape(U,[],nphi);
[I,J] = max(y');
size(J)
size(astruc.points)
p = patch('Faces',astruc.t,'Vertices',astruc.points,'FaceVertexCData',I(:),'FaceColor','interp','EdgeColor','none');

view(2)
axis equal
axis off
toc




