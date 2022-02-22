function result = interferencevectorvectorised(M,Nm,K,cdot,mygamman,mybetan,pilotseq,index)
% The vectorized implementation of interferencevector function
otherusers=1:K~=index;
inteferencevector=mygamman(:,otherusers)./mybetan(:,otherusers).*(repmat(mybetan(:,index),1,K-1));
inteferencevector = cdot(:,index)'*inteferencevector;
inteferencevector=inteferencevector.*abs(pilotseq(:,index)'*pilotseq(:,otherusers));
term2=cdot.*sqrt(mygamman.*repmat(mybetan(:,index),1,K)./Nm);
result=[inteferencevector.';term2(:)];

end

