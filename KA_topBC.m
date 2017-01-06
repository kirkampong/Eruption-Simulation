function residual = KA_topBC( u0,p0,K,H,opt )
% vary u0 until residual=0 (shooting method)
  % to enforce choked flow condition at top boundary
  %
  % residual function is defined piecewise:
  % let u=velocity at top of conduit, 
  % c=sound speed at top of conduit
  %
  % for u<c, residual = u-c<0
  %
  % but if flow becomes sonic below top of conduit,
  % then ODE system (as formulated) becomes singular;
  % solution is only given until some depth |y|, which
  % decreases toward zero as choked flow condition is met
  %
  % so define residual |y|>0 in this case 
  
  q0 = [p0; u0]; % initial conditions
  sol = ode15s(@KA_eruptODE,[-H 0],q0,opt,K); % solve ODE system
  depth = min(-sol.x); % minimum depth at which solution is nonsingular
  if depth>0 % flow becomes sonic below vent, so decrease u0
    residual = depth; % penalize distance from sonic point to vent
  else
    p = sol.y(1,end); % pressure at top of conduit
    u = sol.y(2,end); % velocity at top of conduit
    [~,~,c] = KA_eos(p,K); % sound speed at top of conduit
    residual = u-c; % residual is deviation from choked flow at vent
  end



