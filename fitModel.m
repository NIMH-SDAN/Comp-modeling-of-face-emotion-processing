function KS = fitModel(zr,a,t0,v1,v15,data)

sigma = .1; %root diffusion rate
%zr: starting point relative to a
%a: threshold separation
%t0: nondecision time
%v1,v15 = drift rates for stimuli 1 and 15
%data: [ntrials x 3]
 %data(:,1): stimulus value, in 1:15
 %data(:,2): response, as +/-1
 %data(:,3): RT, in seconds
 %data(:,4): empirical survival probability of each RT within its corresponding stimulus and response values 
 %data(:,5): empirical left-survival probability of each RT within its corresponding stimulus and response values 

G = zeros(size(data,1),1); %model's survivor probability on each trial: probability of observed response at or longer than observed RT

%reparameterize so that upper threshold corresponds to subject's response on each trial
v = (v1 + (v15-v1)*(data(:,1)-1)/14) .* data(:,2); %drift rate on each trial relative to direction of response, [ntrials x 1]
d = a*(1/2-(zr-1/2)*data(:,2)); %distance to winning bound on each trial, [ntrials x 1]

%handle trials with rt < t0
short = data(:,3)<t0; %Boolean, [ntrials x 1]
G(short) = (1-exp(-2*v(short).*(a-d(short))/sigma^2))./(1-exp(-2*v(short)*a/sigma^2)); %model's probability of upper response
G(short & v==0) = 1 - d(short & v==0)/a; %previous expression is udnefined when drift is zero

%now handle long trials (i.e. with rt ≥ t0)
vl = v(~short); %drift on long trials
dl = d(~short); %distance to bound on long trials
t = data(~short,3)-t0; %decision time on long trials

%find maximal number of summands needed for approximation below
delta = 1e-10; %desired accuracy of approximation
A = pi^2*sigma^2*t/(2*a^2); %[ntrials x 1]
kmax=floor(sqrt(max((2*vl.*dl-vl.^2.*t)/(2*sigma^2)-log((1-exp(-A))*pi*delta/2)./A))); %number of summands, scalar

%calculate survivor probabilities on long trials
k = 1:kmax; %summation index, 1xK
B = k.^2*pi^2*sigma^2/(2*a^2) + vl.^2/(2*sigma^2); %[1 x kmax] + [ntrials x 1] = [ntrials x kmax]
C = sin(pi*dl*k/a); %[ntrials x kmax]
D = exp(vl.*dl/sigma^2); %[ntrials x 1]
G(~short) = pi*sigma^2/a^2 * D .* ((C.*exp(-B.*t)./B)*k'); %survivor probabilities on long trials

%plot(data(:,2).*data(:,3),G,'.')

%calculate KS statistic
ks = zeros(15,1); %KS stat for each stimulus
for i=1:15
    stimTrials = data(:,1)==i; %Boolean indicating trials for this stimulus
    ks(i) = max([G(stimTrials)-data(stimTrials,4);data(stimTrials,5)-G(stimTrials)]); %max of over- and under-predictions
end
KS = sum(ks); %sum over stimuli for final fit measure

%for tracking progress of search routine; remove this block when fitting many datasets
global searchlog
searchlog(end+1,:) = [zr,a,t0,v1,v15,KS]; %logging parameters tested by optimization routine, and fit value
if mod(size(searchlog,1),200)==0
    fprintf([num2str(size(searchlog,1)) sprintf(' %f',searchlog(end,:)) '\n']) %display every 200th result
end