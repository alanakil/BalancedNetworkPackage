%% Connectivity of the balanced network %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Recurrent net connection probabilities
P=[0.05 0.05; 0.05 0.05];

% Ffwd connection probs
Px=[.05; .05];

% Mean connection strengths between each cell type pair
Jm=[50 -300; 225 -500]/sqrt(N);
Jxm=[180; 135]/sqrt(N);

% Build mean field matrices
Q=[qe qi; qe qi];
Qf=[qf; qf];
W=P.*(Jm*sqrt(N)).*Q;
Wx=Px.*(Jxm*sqrt(N)).*Qf;

% Generate connectivity matrices
tic
J=[Jm(1,1)*binornd(1,P(1,1),Ne,Ne) Jm(1,2)*binornd(1,P(1,2),Ne,Ni); ...
   Jm(2,1)*binornd(1,P(2,1),Ni,Ne) Jm(2,2)*binornd(1,P(2,2),Ni,Ni)];
Jx=[Jxm(1)*binornd(1,Px(1),Ne,Nx); Jxm(2)*binornd(1,Px(2),Ni,Nx)];


tGen=toc;
disp(sprintf('\nTime to generate connections: %.2f sec',tGen))
