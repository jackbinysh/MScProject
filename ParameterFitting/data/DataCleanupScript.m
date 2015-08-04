clear all;
% a script to clean up the data, and split it into parts
% the cleanup will remove an (assumed) linear drift from some of the data

% this data has a slow drift in it. To clean it up, I'm going to fit an
% exponential decay - we have experimetal evidence that baseline 
% flourescence is 9, so I will subtract this off the data, and fit the
% exponential to the remainder
z0 = 9;

%%% 13_9 %%%
time = load('./InitialExperimentalData/time13_9.csv');
data = load ('./InitialExperimentalData/data13_9.csv');
meandata = mean(data,2);

% find the minima of the oscillations
[x,locs] = findpeaks(-meandata);
%chip away the false minimia
tempmeandata = meandata(locs);
temptime = time(locs);  
while ~isempty(findpeaks(tempmeandata))
    [x,locs] = findpeaks(tempmeandata);
    tempmeandata(locs) = [];
    temptime(locs) = [];  
end

% fit an exponential decay through these
f = fit(temptime,tempmeandata-z0,'exp1');
data = (data-repmat(f.a * exp(f.b*time) ,1,size(data,2)));

% write cleaned data to file
csvwrite('./CleanedData/time13_9.csv',time)
csvwrite('./CleanedData/data13_9.csv',data)

%%% 14_7 %%%
time = load('./InitialExperimentalData/time14_7.csv');
data = load ('./InitialExperimentalData/data14_7.csv');

% this dataset has a few timeseries clearly wrong, I remove these by hand 
data = data(:,~max(data < 8));
meandata = mean(data,2);

% then play the same game,fitting an exponential
% find the minima of the oscillations
[x,locs] = findpeaks(-meandata);
%chip away the false minimia
tempmeandata = meandata(locs);
temptime = time(locs);  
while ~isempty(findpeaks(tempmeandata))
    [x,locs] = findpeaks(tempmeandata);
    tempmeandata(locs) = [];
    temptime(locs) = [];  
end
temptime(end) = []; tempmeandata(end) = []; % last point is spurious.

% fit an exponential
f = fit(temptime,tempmeandata-z0,'exp1');
data = (data-repmat(f.a * exp(f.b*time) ,1,size(data,2)));

% write cleaned data to file
csvwrite('./CleanedData/time14_7.csv',time)
csvwrite('./CleanedData/data14_7.csv',data)

%%% 14_9 %%% 

% get some of the oscillatory bit of 14_9
time = load('./InitialExperimentalData/time14_9.csv');
data = load ('./InitialExperimentalData/data14_9.csv');
time=time(1:90);
data=data(1:90,:);
% write cleaned data to file
csvwrite('./CleanedData/time14_9_1.csv',time)
csvwrite('./CleanedData/data14_9_1.csv',data)

% get the forced but no response portion
time = load('./InitialExperimentalData/time14_9.csv');
data = load ('./InitialExperimentalData/data14_9.csv');
time=time(109:173);
data=data(109:173,:);
csvwrite('./CleanedData/time14_9_2.csv',time)
csvwrite('./CleanedData/data14_9_2.csv',data)
        
 % get the negative control section
time = load('./InitialExperimentalData/time14_9.csv');
data = load ('./InitialExperimentalData/data14_9.csv');
time=time(259:300);
data=data(259:300,:);
csvwrite('./CleanedData/time14_9_3.csv',time)
csvwrite('./CleanedData/data14_9_3.csv',data)






