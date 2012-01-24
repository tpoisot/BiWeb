% Add the code to the matlab path
g = genpath('matlab/');
addpath(g);

% Create a bipartite adjacency matrix that will be used as data
matrix = rand(20) > 0.6;

% Create an object of the class BiWeb. By just creating it,
% modularity and nestedness will be calculated;
w = BiWeb(matrix);

%Accesing the values of nestedness
nodf = w.nestedness.nodf;
nodf_up = w.nestedness.nodf_up;
nodf_low = w.nestedness.nodf_low;

fprintf('Values of Nestedness using NODF algorithm:');
fprintf('\tNODF:                 \t%f', nodf);
fprintf('\tNODF in rows:         \t%f', nod_up);
fprintf('\tNODF in columns:      \t%f', nod_low);

%Accesing the values of modularity using LP&BRIM algorithm

N = w.modules.N;
Qb = w.nestedness.nodf_up;
Qr = w.nestedness.nodf_low;

fprintf('Values of Modularity using LP&BRIM algorithm:');
fprintf('\tNumber of modules:    \t%f', N);
fprintf('\tBipartite modularity: \t%f', Qb);
fprintf('\tRealized modularity:  \t%f', Qr);
