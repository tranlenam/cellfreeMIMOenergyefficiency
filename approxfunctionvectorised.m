function result = approxfunctionvectorised(M,Nm,K,Gammaan,BETAAn,Phii_cf,Pd,etamatrix,u,etamatrix_n,u_n,thetaa,index)
% The vectorized implementation of approxfunction
f1=Pd*(Nm^2)*norm(interferencevectorvectorised(M,Nm,K,etamatrix_n,Gammaan,BETAAn,Phii_cf,index))^2 + 1 +...
    Pd*(Nm^2)*norm(Gammaan(:,index)'*etamatrix_n(:,index))^2;
f=thetaa*f1/u_n(index) - (f1/(u_n(index))^2)*(u(index) - thetaa*u_n(index));

bargamk =(Gammaan./BETAAn.*repmat(BETAAn(:,index),1,K));

t1=sum(etamatrix_n'.*bargamk',2).*abs(Phii_cf'*Phii_cf(:,index));
t2 = sum((etamatrix - thetaa*etamatrix_n)'.*(bargamk)',2).*abs(Phii_cf'*Phii_cf(:,index));

f=f+(2*Pd/u_n(index))*(Nm^2)*sum(t1.*t2);

Dk2 = Gammaan.*repmat(BETAAn(:,index),1,K);

f= f + (2*Pd/u_n(index))*Nm*sum(diag((etamatrix_n.*(Gammaan.*repmat(BETAAn(:,index),1,K)))'...
    *(etamatrix - thetaa*etamatrix_n)));

result=f;

end