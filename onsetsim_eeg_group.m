% onsetsim

corepath = pwd;
addpath([corepath,'/functions'])
ath = 0.05; % arbitrary threshold value for ttest
pth = 1 - ath; % permutation threshold
Nperm = 2000;

%% Make template

% true onset = 160 ms, F=81, max at F=126
true_onset = 160;
Xf = 0:2:500; % 500 Hz
Nf = length(Xf);
temp1 = zeros(1,Nf);
erp = normpdf(linspace(-1.5,1.5,93), 0, 1);
erp = erp - min(erp);
erp = erp / max(erp);

%% ---------------------------------------------------------
%% Simulation: true positives, t^2, Nt=50, var=1

% 20 participants. Group estimate = median of 20 onsets.
% On each iteration, each participant has a random onset of 150-170 ms, drawn from uniform distribution. 

rng(666) % set random seed

Nsim = 10000;
Nperm = 2000;

Nt = 50; % sample size per condition
outvar = 1; % single-trial noise variance
srate = 500; % sampling rate in Hz
Np = 20; % number of participants
ronset = 150:2:170; % random onset for each participant

onset_cluster = zeros(Np, Nsim);
onset_max = zeros(Np, Nsim);
onset_fdr = zeros(Np, Nsim);
cohend = zeros(Np, Nsim);
onset_cp = zeros(Np, Nsim);

% template + noise
cond1 = zeros(Nf, Nt);
cond2 = zeros(Nf, Nt);

for iter = 1:Nsim
    
    if rem(iter,500) == 0
        disp(['onsetsim tp t^2, Nt=50, var=1, EEG noise, group estimation: ',num2str(iter),' / ',num2str(Nsim)])
    end
    
    for P = 1:Np % for each participant
        
        ponset = datasample(ronset, 1); % get random onset
        st = find(Xf==ponset);
        temp2 = [zeros(1,st-2), erp, zeros(1,Nf-st-length(erp)+2)];
%         figure;plot(Xf,temp2);hold on;plot([ponset ponset],[0 1])
                
        for T = 1:Nt
            cond1(:,T) = temp1 + noise(Nf, 1, srate, outvar);
            cond2(:,T) = temp2 + noise(Nf, 1, srate, outvar);
        end
             
        % t-test
        tval = limo_ttest_light(2, cond1, cond2);
        t2 = tval.^2;
        
        % change point:find the two points at which 
        % the mean and SD change the most;
        % onset = first of these two points.
        res = findchangepts(t2, 'Statistic', 'std', 'MaxNumChanges', 2);
        try
            onset_cp(P,iter) = Xf(res(1));
        catch
            onset_cp(P,iter) = NaN;
        end
        
        % get permutation estimates
        t2_perm = zeros(Nf,Nperm);
        erpall = cat(2,cond1,cond2);
        
        for perm_iter = 1:Nperm
            
            perm_trials = randperm(Nt + Nt);
            perm_cond1 = erpall(:,perm_trials(1:Nt));
            perm_cond2 = erpall(:,perm_trials(Nt+1:Nt+Nt));
            tval = limo_ttest_light(2,perm_cond1,perm_cond2);
            t2_perm(:,perm_iter) = tval.^2;
            
        end
        
        maxt2_perm = max(t2_perm, [], 1);
        pval_perm = zeros(Nf,1);
        for F = 1:Nf
            pval_perm(F) = (sum(t2_perm(F,:) > t2(F)) + 1) / (Nperm + 1);
        end
        
        % CLUSTER onset estimate
        
        % get threshold
        perm_th = prctile(t2_perm, pth*100, 2); % univariate thresholds
        % because there is no p value for max t^2,
        % we use the univariate permutation distributions to threshold themselves
        % fake p values: 0 if above threshold, 1 otherwise
        pval = t2_perm <= repmat(perm_th, 1, Nperm);
        % threshold permutation distribution
        tmp = t2_perm;
        th = limo_ecluster_make(tmp, pval, ath);
        % threshold T2 results
        sigcluster = limo_ecluster_test(t2, t2 < perm_th, th, ath);
        % find onset
        try
            onset_cluster(P,iter) = find_onset(sigcluster.elec_mask, Xf, 1);
        catch
            onset_cluster(P,iter) = NaN;
        end
        
        % MAX onset estimate
        max_perm_th = prctile(maxt2_perm, pth*100); % MAX stat threshold
        try
            onset_max(P,iter) = find_onset(t2 > max_perm_th, Xf, 1);
        catch
            onset_max(P,iter) = NaN;
        end
        
        % FDR onset estimate
        [pID,pN] = FDR(pval_perm, ath);
        try
            onset_fdr(P,iter) = find_onset(pval_perm < pID, Xf, 1);
        catch
            onset_fdr(P,iter) = NaN;
        end
        
    end
    
end

save([corepath,'/onsetsim_t2_n50_var1_eegnoise_group'], ...
    'onset_cluster', 'onset_max', 'onset_fdr', 'cohend', ...
    'onset_cp', 'Nsim', 'Nperm')

%% Simulation results: true positives, t^2 -- one gamma/outvar/Nt
% Check results for one combination of variables.

Nsim = 10000;

true_onset = 160;
load([corepath,'/data/onsetsim_t2_n50_var1_eegnoise_group'])

% median onsets across participants 
oc = median(onset_cluster, 1, 'omitnan');
om = median(onset_max, 1, 'omitnan');
of = median(onset_fdr, 1, 'omitnan');
ocp = median(onset_cp, 1, 'omitnan');
mdo = [median(oc) median(om) median(of) median(ocp)];

disp('-----------------------------------------------------')
disp('Median onsets:')
disp(['cluster = ',num2str(round(median(oc))),' ms'])
disp(['max = ',num2str(round(median(om))),' ms'])
disp(['fdr = ',num2str(round(median(of))),' ms'])
disp(['cp = ',num2str(round(median(ocp))),' ms'])
disp('-----------------------------------------------------')
disp('Mean absolute error (MAE):')
disp(['cluster = ',num2str(round(mean(abs(oc-true_onset)))),' ms'])
disp(['max = ',num2str(round(mean(abs(om-true_onset)))),' ms'])
disp(['fdr = ',num2str(round(mean(abs(of-true_onset)))),' ms'])
disp(['cp = ',num2str(round(mean(abs(ocp-true_onset)))),' ms'])
disp('-----------------------------------------------------')
disp('Bias:')
disp(['cluster = ',num2str(round(median(oc)-true_onset)),' ms'])
disp(['max = ',num2str(round(median(om)-true_onset)),' ms'])
disp(['fdr = ',num2str(round(median(of)-true_onset)),' ms'])
disp(['cp = ',num2str(round(median(ocp)-true_onset)),' ms'])
disp('-----------------------------------------------------')
disp('Underestimations of at least 40 ms:')
disp(['cluster = ',num2str(round(100*mean((oc-true_onset)<= -40),1)),' %'])
disp(['max = ',num2str(round(100*mean((om-true_onset)<= -40),1)),' %'])
disp(['fdr = ',num2str(round(100*mean((of-true_onset)<= -40),1)),' %'])
disp(['cp = ',num2str(round(100*mean((ocp-true_onset)<= -40),1)),' %'])
disp('-----------------------------------------------------')
disp('Proportion too early:')
disp(['cluster = ',num2str(round(100*mean((oc-true_onset)< 0),1)),' %'])
disp(['max = ',num2str(round(100*mean((om-true_onset)< 0),1)),' %'])
disp(['fdr = ',num2str(round(100*mean((of-true_onset)< 0),1)),' %'])
disp(['cp = ',num2str(round(100*mean((ocp-true_onset)< 0),1)),' %'])
disp('-----------------------------------------------------')
disp('Variance:')
disp(['cluster = ',num2str(round(var(oc,1))),' ms'])
disp(['max = ',num2str(round(var(om,1))),' ms'])
disp(['fdr = ',num2str(round(var(of,1))),' ms'])
disp(['cp = ',num2str(round(var(ocp,1))),' ms'])

figure('NumberTitle','off', 'Name', ['nt=',num2str(Nt),', var=',num2str(outvar),', EEG noise'])

subplot(1,4,1)
[f,x] = ksdensity(oc);
plot(x,f,'k')
% histogram(oc, 'Normalization', 'count')
title(['CLUSTER median onset = ',num2str(round(mdo_c))])

subplot(1,4,2)
[f,x] = ksdensity(om);
plot(x,f,'k')
% histogram(om, 'Normalization', 'count')
title(['MAX median onset = ',num2str(round(mdo_m))])

subplot(1,4,3)
[f,x] = ksdensity(of);
plot(x,f,'k')
% histogram(of, 'Normalization', 'count')
title(['FDR median onset = ',num2str(round(mdo_f))])

subplot(1,4,4)
[f,x] = ksdensity(ocp);
plot(x,f,'k')
title(['CP median onset = ',num2str(round(mdo_cp))])

for sub = 1:4
    subplot(1,4,sub)
    hold on
    xlim([0 500])
    v = axis;
    plot([150 150], [v(3) v(4)], 'k', 'LineWidth', 1)
    plot([mdo(sub) mdo(sub)], [v(3) v(4)], 'k--', 'LineWidth', 1)
end




