function [c_n,u_n,t_n] = generateinitialpoint(M,K,N,mygamma,mybeta,rho_d,pilotseq,RateQoS)
% This function is used to generate a feasible initial point to start
% Algorithm. The procedure is as follows
% 1. Generate feasible solution to (24) by solving a feasibility SOCP
% problem
% 2. If a feasible solution is found, then calculate the in
opts=sdpsettings('solver','mosek','verbose',0,'dualize',0); % set internal solver to mosek
c = sdpvar(M,K,'full') ;
F = [c(:)>=0]; % (24d)
F = [F,cone([1/sqrt(N)*ones(1,M);(sqrt(mygamma).*c)'])]; %(24c),
for iUser = 1:K
    F=[F,cone([(1/sqrt(2^( RateQoS(iUser)) - 1))*c(:,iUser)'*(sqrt(rho_d)*mygamma(:,iUser));...
        %            interference_sca_vectorised(M,N,K,cdot,sqrt(rho_d)*mygamma,sqrt(rho_d)*mybeta,pilotseq,iUser);mytheta/N])]; %(33d)
        interferencevector(M,N,K,c,sqrt(rho_d)*mygamma,sqrt(rho_d)*mybeta,pilotseq,iUser);1/N])]; %(33d)
end
diagnotics = optimize(F,[],opts); % solve the feasibility problem
if diagnotics.problem==0 % the problem is feasible
    c_n = value(c);
    t_n = RateQoS;
    u_n = 2.^t_n;
else % problem is infeasible, assign NaN to all variables
    t_n = NaN;
    u_n = NaN;
    c_n = NaN;
end
end

