function SensitivityMatrix = SensitivityMatrix(theta,Dataset)
	
	% read in the dataset, and get the integration times
    Times = csvread(strcat('./data/CleanedData/time',Dataset,'.csv'));

	% choose the model to use 
	cd('./model1')
	
	% options for solver
	x0 = InitialState(theta,0,1);
	options = odeset('Jacobian',@Jacobian,'MaxStep',10);
	
	for i = 1:length(theta);
	
		% amount to perturb by
		dtheta = 0.01*theta(i);
		
		% get the upper and lower perturbed thetas
		Uppertheta = theta; Lowertheta = theta;
		Uppertheta(i) = theta(i) + dtheta;
		Lowertheta(i) = theta(i) - dtheta;
	
		% simulate
		[T,UpperPrediction] = ode15s(@RibodynamicsModel,[0:0.1:max(Times)],x0,options,Uppertheta,Dataset);
		[T,LowerPrediction] = ode15s(@RibodynamicsModel,[0:0.1:max(Times)],x0,options,Lowertheta,Dataset);
		
		% compute numerical derivatives
		SensitivityMatrix(:,i) =  ( UpperPrediction(:,6) - LowerPrediction(:,6) )/ (2*dtheta);
	end
	
	cd('..')
end
