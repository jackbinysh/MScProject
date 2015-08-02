
format long g;

% decide which dataset we are using
Dataset = '16_8';

% decide which run we use: if we leave it empty, an average is taken

% Initial Guess and bounds. Thesse  are currently set to the best guess
% from GArun_30_07_2015.
% list of variables {'f_srna';'k_on';'k_off';'k_hyb';'delta_m';'delta_s';'delta_c';'mu';'beta';'c'};
InitialthetaLB = [1;1000;100000;1;1;1;1;0.00001;0.0001;1];
% from GArun_30_07_2015.
Initialtheta = [ 51.6568541032315
          333746.163658966
          35923992.3356378
          78.0092377634537
          52.5614560200565
          40.2395554589868
          70.8626480270463
        0.0179550155596016
         0.184074335687385
          401.122143648984 ];

InitialthetaUB = [1e3;1e6;1e8;1000;1000;1000;1000;1;10;10000];

% plug into the wrapper function for the fitter
[xmin,fmin,counteval,stopflag,out] = ParameterFit(Initialtheta, InitialthetaLB, InitialthetaUB, Dataset,RunNumber);
