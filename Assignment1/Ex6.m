L=7; % number of cells
Nb=2^L-1; % m-sequence length
pnSequence = comm.PNSequence('Polynomial',[7 3 0], ...
'SamplesPerFrame',Nb,'InitialConditions',[1 0 1 0 1 0 1]);
x1 = pnSequence()'; % binary m-sequence