function modelEstimation(inFile, outFile)

%Estimate parameters of diffusion model on fIBT data
%assuming drift as a linear function of stimulus

%model parameters
 %zr: starting point relative to a
 %a: threshold separation
 %t0: nondecision time
 %v1,v15 = drift rates for stimuli 1 and 15
 
%data structure
 %data(:,1): stimulus value, in 1:15
 %data(:,2): response, as 0/1 where 1 is the nominal response for stimulus 15
 %data(:,3): RT, in seconds
 %code currently assumes minimum of 1 trial for each stimulus
 
[data] =csvread(inFile);

global searchlog
searchlog = [];

%count presentations of each stimulus
trialsPerStim = zeros(15,1); 
for s=1:15
    trialsPerStim(s) = sum(data(:,1)==s);
end

%restructure data
dataReformat = [data(:,1),2*data(:,2)-1,data(:,3)]; %recode response as +/-1
for i=1:size(data,1) %loop trials
    rank = sum(data(:,1)==data(i,1) & data(:,2)==data(i,2) & data(:,3)>=data(i,3)); %reverse rank of current RT among all trials with same stimulus and response (i.e. transformed cumulative frequency)
    dataReformat(i,4) = (rank-1)/trialsPerStim(data(i,1)); %empirical (right-)survivor probability
    dataReformat(i,5) = rank/trialsPerStim(data(i,1)); %empirical left-survivor probability
end

%estimate model parameters, coded as theta = [zr,a,t0,v1,v15]
f = @(theta) fitModel(theta(1),theta(2),theta(3),theta(4),theta(5),dataReformat); %objective function
LB = [.1,.01,.1,-2,-2]; %lower bounds for parameters
UB = [.9,  1, 1, 2, 2]; %upper bounds for parameters
options = optimoptions('ga','MaxGenerations',5000);
[theta,KS] = ga(f,5,[],[],[],[],LB,UB,[],[],options); %optimize parameters

%extract parameters
zr = theta(1); %response bias; zr>.5 is bias toward response 1 vs response 0
a = theta(2); %threshold separation
t0 = theta(3); %nondecision time
sens = (theta(5)-theta(4))/14; %sensitivity: change in drift per unit change in stimulus
indiff = 1 - theta(4)/sens; %indifference point: stimulus value producing zero drift

save(outFile,'theta','KS','zr','a','t0','sens','indiff');

end
