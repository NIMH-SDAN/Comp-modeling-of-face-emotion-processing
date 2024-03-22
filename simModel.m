function data = simModel(zr,a,t0,v1,v15)
%simulate diffusion model to create synthetic dataset
%zr = starting point relative to threshold separation
%a = upper threshold, 0 = lower threshold
%t0 = nondecision time
%v1,v15 = drift rates for stimuli 1 and 15

sigma = .1; %root diffusion rate
delta = .0001; %temporal resolution for simulating diffusion process
nstim = 15; %stimulus values range 1:nstim
trialsPerStim = 1000; %each stimulus presented trialsPerStim times

data = zeros(nstim*trialsPerStim,3); %[stim response rt] on each trial; response as 0/1, rt in seconds
for s=1:nstim %loop stimuli
    v = v1 + (v15-v1)*(s-1)/14; %drift rate for this stimulus
    for t=1:trialsPerStim %loop trials for this stimulus
        e = zr*a; %initialize evidence at starting point
        rt = t0; %initialize RT at nondecision time
        while 1 %run diffusion process until a threshold is reached
            rt = rt+delta; %increment time
            e = e + v*delta + randn(1)*sigma*sqrt(delta); %increment evidence
            if e<0 %lower threshold crossed
                data((s-1)*trialsPerStim+t,:) = [s,0,rt]; %record data for this trial
                break %end trial
            elseif e>a %upper threshold crossed
                data((s-1)*trialsPerStim+t,:) = [s,1,rt]; %record data for this trial
                break %end trial
            end
        end
    end
end
