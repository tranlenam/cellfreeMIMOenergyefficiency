function result = approxfunction(M,N,K,mygamma,mybeta,pilotseq,rho_d,cdot,udot,c_n,u_n,mytheta,iUser)
% This function returns the approximate function defined in (37)

f1 = rho_d*(N^2)*norm(interferencevector(M,N,K,c_n,mygamma,mybeta,pilotseq,iUser))^2 ...
    + 1 + rho_d*(N^2)*norm(mygamma(:,iUser)'*c_n(:,iUser))^2; % refer to (30)

f= mytheta*f1/u_n(iUser) - (f1/(u_n(iUser))^2)*(udot(iUser) - mytheta*u_n(iUser));

for k=1:K
    mygammabar = ((mygamma(:,k)./mybeta(:,k)).*mybeta(:,iUser))*abs(pilotseq(:,k)'*pilotseq(:,iUser));
    Dk2 = diag( mygamma(:,k).*mybeta(:,iUser) );
    f = f + (2*rho_d/u_n(iUser))*c_n(:,k)'*((N^2)*(mygammabar*mygammabar') + N*Dk2)*(cdot(:,k) - mytheta*c_n(:,k));
end

result=f;

end