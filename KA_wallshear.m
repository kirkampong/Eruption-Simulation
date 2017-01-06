function tau = KA_wallshear(rho,phi,rhoc,rhol,rhoe,rhod,Xd,Xe,u,K,c,beta)

switch K.opt1
    case 'dynamic'
        % effect of volatiles on viscosity [Hess & Dingwell 1996]
        logmu = (-3.545 + 0.833*log(Xd*100)) + ...
            (9601 - 2368*log(Xd*100))./((K.T-(195.7 + 32.25*log(Xd*100))));
        mu  = 10.^logmu;
        
    case 'constant'
        mu = 1e6; % constant viscosity
        
    otherwise
        error('Invalid volatile viscosity effect option entered');
        
end

switch K.opt2
    case 'standard'
        % Effect of crystals on viscosity [Costa 2005]
        phi_c = ((K.PHI./rhoc)./(1-K.n0-K.PHI))./...  %crystal volume fraction
            ((Xe./rhoe)+(Xd./rhod)+(1./rhol)+((K.PHI./rhoc)./(1-K.n0-K.PHI)));
        
        a = 0.9995; B = 0.4; Y = 1;   % parameters
        x = 0.5.*sqrt(pi).*phi_c.*(1+(B./((1-phi_c).^Y)));
        
        mu_rel = (1-a.*erf(x)).^(-2.5/a);  % relative viscosity
        mu = mu.*mu_rel; %Effective viscosity
        
    case 'improved'
        % Effect of crystals on viscosity [Improved Costa 2005]
        phi_c = ((K.PHI./rhoc)./(1-K.n0-K.PHI))./...  %crystal volume fraction
            ((Xe./rhoe)+(Xd./rhod)+(1./rhol)+((K.PHI./rhoc)./(1-K.n0-K.PHI)));
        
        a = 0.999916 ; B = 2.5; s = 16.9386; Y = 3.98937; phi_m = 0.673;
        x = ((sqrt(pi).*phi_c)/(2*a.*phi_m)).*(1+(phi_c./phi_m).^Y);
        
        mu_rel = (1+(phi_c./phi_m).^s)./((1-a.*erf(x)).^(B.*phi_m));% relative viscosity
        mu = mu.*mu_rel; %Effective viscosity
        
    otherwise
        ; % no enhancement due to crystals, or...
        % error('Invalid crystal viscosity effect option entered'); % to throw error
        
end

switch K.opt3
    case 'shear'
        %shear fragmentation model
        
        %strain rate below fragmentation depth:
        dudy1 = -u.*beta.*-((K.f0.*rho.*u.^2+8.*mu.*u./K.r)./K.r+rho.*K.g)./(1-(u./c).^2);
        
        %strain rate above fragmentation depth:
        dudy2 = -u.*beta.*-((K.f0.*rho.*u.^2)./K.r+rho.*K.g)./(1-(u./c).^2);
        
        if (dudy1 > (0.01*(K.Kl./mu)))
            tau = K.f0.*0.5.*rho.*u.^2;
        else
            tau = K.f0.*0.5.*rho.*u.^2 + 4.*mu.*u./K.r; %laminar
        end
        
    case 'gas_fraction'
        %critical gas fraction fragmentation model
        tau = K.f0.*0.5*rho.*u.^2; % turbulent wall shear stress
        % added laminar contribution below fragementation depth
        tau(phi<K.phi0) = tau(phi<K.phi0) + 4*mu(phi<K.phi0).*u(phi<K.phi0)./K.r;
        
    otherwise
        error('Invalid fragemntation option entered');
        
end

end


