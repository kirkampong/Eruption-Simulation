function [rho,phi,c,beta,rhoc,rhol,rhoe,rhod,Xd,Xe,ne,nd] = KA_eos(p,K )

    
  p0 = (K.X0/K.s)^(1/K.m); % p>p0 => all gas in solution
  
  % solubility law
  Xd = min(K.s*p.^K.m,K.X0); % dissolved gas mass concentration
  Xe = K.X0 - Xd; % exsolved gas mass concentration 
  
  %ne, nd
  ne = Xe./(1+Xe);
  nd = Xd./(1+Xd);
  
  % liquid
  rhol = K.rhol0*(1+(p-p0)/K.Kl); % liquid density
  drholdp = K.rhol0/K.Kl;
   
  % gas   
  rhoe = p/(K.R*K.T); % exsolved gas density
  drhoedp = 1/(K.R*K.T);
  rhod = K.rhol0*(1+(p-p0)/K.Kl); % dissolved gas density
  drhoddp = K.rhol0/K.Kl;
  
  dXddp = K.m*K.s.*p.^(K.m-1); dXddp(ne==0) = 0;
  dXedp = -dXddp; dXedp(ne==0) = 0; 
  
  % crystals
  rhoc = K.rhoc0*(1+(p-p0)/K.Kl_c); % crystal density
  drhocdp = K.rhoc0/K.Kl_c;
  
  dPHIdp = 0; 
  
  % mixture density
  rho_inverse = (((1-K.PHI)./(1+K.X0)).*((1./rhol)+(Xd./rhod)+(Xe./rhoe)))+...
                (K.PHI./rhoc);
  rho = 1./rho_inverse;
  
  % mixture compressibility
  % beta = (1/rho)*(drho/dp) = -rho*(d(RHO)/dp) where RHO = rho_inverse
  dRHOdrhol = -(1-K.PHI)./((1+K.X0).*rhol.^2);  
  dRHOdrhoe = -Xe.*(1-K.PHI)./((1+K.X0).*rhoe.^2);
  dRHOdrhod = -Xd.*(1-K.PHI)./((1+K.X0).*rhod.^2);
  dRHOdrhoc = -K.PHI./rhoc.^2;
  dRHOdXe = (1-K.PHI)./((1+K.X0).*rhoe);  
  dRHOdXd = (1-K.PHI)./((1+K.X0).*rhod); 
  dRHOdPHI = -((1./rhol)+(Xe./rhoe)+(Xd./rhod))./(1+K.X0) + (1./rhoc); 
  
  dRHOdp = (dRHOdrhol.*drholdp)+(dRHOdrhoe.*drhoedp)+(dRHOdrhod.*drhoddp)...
          +(dRHOdrhoc.*drhocdp)+(dRHOdXe.*dXedp)+(dRHOdXd.*dXddp)+...
           (dRHOdPHI.*dPHIdp);
       
  beta = -rho.*dRHOdp;
  
  % sound speed
  c = 1./sqrt(rho.*beta);
  
  % exsolved gas volume fraction
  phi = ((Xe./rhoe))./((Xe./rhoe)+(Xd./rhod)+(1./rhol)...
          +((K.PHI./rhoc)./(1-K.n0-K.PHI)));
  
end

