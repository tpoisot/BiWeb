% Add the code to the matlab path
g = genpath('matlab/');
addpath(g);

% Create an object of the class BiWeb. By just creating it,
% nestedness will be calculated;
w = Reading.CREATE_FROM_MATRIX_WITH_LABELS('sampledata/rodents.web');

%Accesing the values of nestedness
nodf = w.nestedness.nodf;
nodf_up = w.nestedness.nodf_rows;
nodf_low = w.nestedness.nodf_cols;

fprintf('\nValues of Nestedness using NODF algorithm:\n');
fprintf('\tNODF:                 \t%f\n', nodf);
fprintf('\tNODF in rows:         \t%f\n', nodf_up);
fprintf('\tNODF in columns:      \t%f\n', nodf_low);

%Accesing the values of modularity using LP&BRIM algorithm
w.modules.Detect();
N = w.modules.N;
Qb = w.modules.Qb;
Qr = w.modules.Qr;

fprintf('\nValues of Modularity using LP&BRIM algorithm:\n');
fprintf('\tNumber of modules:    \t%f\n', N);
fprintf('\tBipartite modularity: \t%f\n', Qb);
fprintf('\tRealized modularity:  \t%f\n', Qr);


%Print species level properties;
%The properties will be printed to screen and to file rodents-sp
w.printer.SpeciesLevel()

%Create statistics for modularity and nestedness;
n_trials = 10; %Increase this value for a best accurate result
w.tests.DoNulls(@NullModels.NULL_1,n_trials); %Using bernoulli null model
w.tests.Nestedness();
w.tests.Modularity();
w.printer.NetworkLevel(); %Add a row with the statistics to networklevel file
w.tests.DoNulls(@NullModels.NULL_2,n_trials); %Using average in columns+rows
w.tests.Nestedness();
w.tests.Modularity();
w.printer.NetworkLevel();
w.tests.DoNulls(@NullModels.NULL_3COL,n_trials); %Using average in columns
w.tests.Nestedness();
w.tests.Modularity();
w.printer.NetworkLevel();
w.tests.DoNulls(@NullModels.NULL_3ROW,n_trials); %Using average in rows
w.tests.Nestedness();
w.tests.Modularity();
w.printer.NetworkLevel();

%Plot network as matrices
figure(1);
w.plotter.PlotMatrix();
figure(2);
w.plotter.PlotNestedMatrix();
figure(3);
w.plotter.PlotModularMatrix();

%Plot network as beads (graph)
figure(4);
w.plotter.PlotBMatrix();
figure(5);
w.plotter.PlotBNestedMatrix();
figure(6);
w.plotter.PlotBModularMatrix();
