function BETAA = getslowfading(channelparams)
%{
Consider a rectangular area with DxD m^2
M distributed APs serves K terminals, they all randomly located in the
area
%}

M = channelparams.nAPs;
K = channelparams.nUsers;
L = channelparams.pathloss;
D = channelparams.dim;
sigma_shd = channelparams.shadowdev;
d0 = channelparams.refdist0;
d1 = channelparams.refdist1;


%Randomly locations of M APs:
AP=unifrnd(-D/2,D/2,M,2);
%Wrapped around (8 neighbor cells)
D1=zeros(M,2);
D1(:,1)=D1(:,1)+ D*ones(M,1);
AP1=AP+D1;

D2=zeros(M,2);
D2(:,2)=D2(:,2)+ D*ones(M,1);
AP2=AP+D2;

D3=zeros(M,2);
D3(:,1)=D3(:,1)- D*ones(M,1);
AP3=AP+D3;

D4=zeros(M,2);
D4(:,2)=D4(:,2)- D*ones(M,1);
AP4=AP+D4;

D5=zeros(M,2);
D5(:,1)=D5(:,1)+ D*ones(M,1);
D5(:,2)=D5(:,2)- D*ones(M,1);
AP5=AP+D5;

D6=zeros(M,2);
D6(:,1)=D6(:,1)- D*ones(M,1);
D6(:,2)=D6(:,2)+ D*ones(M,1);
AP6=AP+D6;

D7=zeros(M,2);
D7=D7+ D*ones(M,2);
AP7=AP+D7;

D8=zeros(M,2);
D8=D8- D*ones(M,2);
AP8=AP+D8;

%Randomly locations of K terminals:
Ter=unifrnd(-D/2,D/2,K,2);



%Create an MxK large-scale coefficients beta_mk
BETAA = zeros(M,K);
dist=zeros(M,K);
for m=1:M
    for k=1:K
        dist(m,k) = min([norm(AP(m,:)-Ter(k,:)), norm(AP1(m,:)-Ter(k,:)),norm(AP2(m,:)-Ter(k,:)),norm(AP3(m,:)-Ter(k,:)),...
            norm(AP4(m,:)-Ter(k,:)),norm(AP5(m,:)-Ter(k,:)),norm(AP6(m,:)-Ter(k,:)),norm(AP7(m,:)-Ter(k,:)),norm(AP8(m,:)-Ter(k,:)) ]); %distance between Terminal k and AP m
        
        if dist(m,k)<d0
            betadB=-L - 35*log10(d1) + 20*log10(d1) - 20*log10(d0);
        elseif ((dist(m,k)>=d0) && (dist(m,k)<=d1))
            betadB= -L - 35*log10(d1) + 20*log10(d1) - 20*log10(dist(m,k));
        else
            betadB = -L - 35*log10(dist(m,k)) + sigma_shd*randn(1,1); %large-scale in dB
        end
        
        BETAA(m,k)=10^(betadB/10);
    end
end

end