function dqdy = KA_eruptODE(y,q,K)

  p = q(1,:); u = q(2,:);
  N = size(q,2); dqdy = zeros(2,N);
  
  [rho,phi,c,beta,rhoc,rhol,rhoe,rhod,Xd,Xe] = KA_eos(p,K);
  
  tau = KA_wallshear(rho,phi,rhoc,rhol,rhoe,rhod,Xd,Xe,u,K,c,beta);
  
  rhs = -(2*tau./K.r+rho*K.g)./(1-(u./c).^2);
  
  dqdy(1,:) = rhs;
  dqdy(2,:) = -u.*beta*rhs;
  
  
end

