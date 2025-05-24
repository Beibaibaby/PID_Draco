function fit_pair_label_allSessions(data_master, srcField, dstField, ...
                                    labelType, rootDir, tauVal)

pairTag = sprintf('%s_%s_%s',srcField,dstField,labelType);
matDir  = fullfile(rootDir,pairTag,'session_mat');  if ~exist(matDir,'dir'), mkdir(matDir); end
figDir  = fullfile(rootDir,pairTag,'session_fig');  if ~exist(figDir,'dir'), mkdir(figDir); end

% ----- analysis parameters (unchanged from your script) -----------------
alignEvt     = 'Align_to_cat_stim_on';
fs           = 40000;
win          = 0.050;
step         = 0.010;
t0           = 0.0;   t1 = 0.3;
centres      = t0:step:t1;
nW           = numel(centres);
subBin       = 0.005;
minTrials    = 30;

optsFIT = struct('tauVal',tauVal, ...        % <<<<<<<<<<<<<<<<<<<<<<
                 'redundancy_measure','I_min', ...
                 'bin_method',{{'eqpop','eqpop','none'}}, ...
                 'n_bins',{{3,3}}, ...
                 'bias','plugin','supressWarnings',true, ...
                 'pid_constrained',true);


% ----- sessions ----------------------------------------------------------
sessIDs = [data_master.Bhv.session_id]';
sessions = unique(sessIDs);     nSess = numel(sessions);

Aunits   = data_master.Neuro.(srcField);
Bunits   = data_master.Neuro.(dstField);
datesA   = parseDates({Aunits.NeuronID});
datesB   = parseDates({Bunits.NeuronID});

for s = 1:nSess
    sd = sessions(s);
    mask = sessIDs==sd;
    T    = data_master.Bhv(mask).Trial_info;
    T    = T(~isnan([T.category]));
    if numel(T)<minTrials, continue; end

    dirLab = [T.direction]';  catLab = [T.category]';  nTr = numel(T);

    A = Aunits(datesA==sd);  B = Bunits(datesB==sd);
    if isempty(A)||isempty(B), continue; end

    dirAB = nan(nW,1); dirBA = dirAB; catAB = dirAB; catBA = dirAB;

    parfor w = 1:nW
        tc = centres(w);
        ts = max(tc-win/2,0); te = ts+win;

        [~,XA,YB] = extract_spike_matrices(T,A,B,fs,[ts,te],alignEvt,subBin);
        if isempty(XA)||isempty(YB),  continue; end
        bins=min(size(XA,2),size(YB,2));
        XA=sum(XA(:,1:bins,:),3); YB=sum(YB(:,1:bins,:),3);

        X1=reshape(XA,[1,bins,nTr]);  Y1=reshape(YB,[1,bins,nTr]);
        opts = optsFIT; opts.tpres={bins}; opts.tau={min(tauVal,bins-1)}; 

        if strcmp(labelType,'direction')
            lab=reshape(dirLab,[1,1,nTr]);
            v=FIT({X1,Y1,lab},{'FIT(A->B;S)','FIT(B->A;S)'},opts);
            dirAB(w)=v{1}; dirBA(w)=v{2};
        else
            lab=reshape(catLab,[1,1,nTr]);
            v=FIT({X1,Y1,lab},{'FIT(A->B;S)','FIT(B->A;S)'},opts);
            catAB(w)=v{1}; catBA(w)=v{2};
        end
    end % parfor

    % save per-session
    sess.time_centers=centres;
    sess.dirAtoB=dirAB; sess.dirBtoA=dirBA;
    sess.catAtoB=catAB; sess.catBtoA=catBA;
    save(fullfile(matDir,sprintf('sess_%d.mat',sd)),'sess');

    % quick png
    tt=centres*1e3;
    if strcmp(labelType,'direction')
        quickPlot(tt,dirAB,dirBA,[pairTag ' DIR'],fullfile(figDir,sprintf('%d_dir.png',sd)));
    else
        quickPlot(tt,catAB,catBA,[pairTag ' CAT'],fullfile(figDir,sprintf('%d_cat.png',sd)));
    end
end
end

% -- helpers --------------------------------------------------------------
function dates=parseDates(ids)
dates=nan(size(ids));
for k=1:numel(ids)
  t=regexp(ids{k},'^(\d{8})_','tokens','once');
  if ~isempty(t), dates(k)=str2double(t{1}); end
end
end
function quickPlot(t,x,y,ttl,fname)
fig=figure('visible','off'); hold on
plot(t,x,'-r','LineWidth',2); plot(t,y,'-b','LineWidth',2);
refline(0,0); xlabel('ms'); ylabel('FIT'); title(ttl); grid on
legend('A→B','B→A'); saveas(fig,fname); close(fig);
end
