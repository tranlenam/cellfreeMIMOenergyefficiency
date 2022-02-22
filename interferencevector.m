function result = interferencevector(M,Nm,K,cdot,mygamman,mybetan,pilotseq,iUser)
% compute a vector that contains the elements relating to cdot in (36d)
temp = [];
for jUser=1:K
    if jUser~=iUser
        bargamj = ((mygamman(:,jUser)./mybetan(:,jUser)).*mybetan(:,iUser))*abs(pilotseq(:,jUser)'*pilotseq(:,iUser));
        temp = [temp,bargamj'*cdot(:,jUser)];
    end
end

for k=1:K
    for m=1:M
        temp = [temp,(1/sqrt(Nm))*cdot(m,k)*sqrt(mygamman(m,k)*mybetan(m,iUser))];
    end
end

result=temp.'; % convert to a column vector
end