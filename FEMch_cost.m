function [valf,gradt] = cost_spherep3(U,opts)

if isfield(opts,'Prob')==0
   Prob = opts;
else
   Prob = opts.Prob;
end

ep   = Prob.ep;
nphi = Prob.nphi;
npt  = Prob.npt;
K    = Prob.K;
M    = Prob.M;
dim  = 2;
pow  = 2*dim/(dim-1);
alpha = Prob.alpha;
fact = Prob.fact;
pn   = Prob.pn;
Ibord = Prob.Ibord;

U = reshape(U,npt,nphi);
v  = U.*(1-U);
ovec = ones(npt,1);
gradt = zeros(size(U));
gradp = gradt;
gvol  = zeros(size(U));
vals  = zeros(nphi,1);

for i=1:nphi
  pers(i) = ep*U(:,i)'*K*U(:,i)+1/ep*v(:,i)'*M*v(:,i);
  gradp = 2*ep*K*U(:,i)+1/ep*2*M*v(:,i).*(1-2*U(:,i));

  vol(i)      =     ovec'*M*U(:,i).^pow;
  gvol(:,i)   = pow*(ovec'*M)'.*U(:,i).^(pow-1);   

  vals(i)                          = pers(i)/(vol(i))^alpha;
  gradt(:,i)                       = (vol(i)^alpha*gradp-alpha*vol(i)^(alpha-1)*gvol(:,i)*pers(i))/vol(i)^(2*alpha);
end

pen = 0;
gradpen = zeros(size(U));


for i=1:nphi
 for j=i+1:nphi
    pen = pen+U(:,i)'.^2*M*U(:,j).^2;
    gradpen(:,i) = gradpen(:,i)+2*U(:,i).*(M*U(:,j).^2);
    gradpen(:,j) = gradpen(:,j)+2*U(:,j).*(M*U(:,i).^2);
 end
end

val = vals;
valf                                = (sum(val.^pn))^(1/pn);
for i=1:nphi
  gradt(:,i)       = gradt(:,i)*(sum(val.^pn))^(1/pn-1)*val(i)^(pn-1);
end

valf  = valf+fact*pen;
gradt = gradt+fact*gradpen;

gradt(Ibord,:) = 0;

gradt                              = gradt(:);





