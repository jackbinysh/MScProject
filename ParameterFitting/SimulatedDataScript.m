% a script to generate some artificial data, which 
% we can test our algorithm on

% artificial parameters
theta = [100;560000;29126800;77;50;0.08;62;0.003;0.003;422];

%artificial initial state
x0 = [120.4994;1000;234;751;0.27;10];

% generate model prediction
Times = csvread('../data/InitialExperimentalData/time16_8.csv');
%MEX
cd('./MexRibodynamics')
ModelWithParams = @(t,x) RibodynamicsModel(t,x', theta);
% get the model prediction
[T,Prediction] = ode23s(ModelWithParams,Times,x0);
cd('..')
Prediction = Prediction(:,6);
% now fuzz the time series with gaussian noise.
%First set a scale for the noise
Scale = 0.1* (max(Prediction(:,6)) - mean(Prediction(:,6)));

% now generate noisy paths
NumSeries = 145;
Noise = normrnd(0,repmat(Scale,length(Prediction),NumSeries));
SimulatedData = repmat(Prediction,1,NumSeries) + Noise;

