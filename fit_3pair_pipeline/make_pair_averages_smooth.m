function make_pair_averages_smooth(resultsRoot, win, method, pct)
% make_pair_averages_smooth  – smoothed mean ±STD, unified y-scale
%
%   resultsRoot : where pair folders live           (default 'results')
%   win         : window length for smoothdata      (default 3 samples)
%   method      : smoothdata method ('movmean' …)   (default 'movmean')
%   pct         : central percentile kept for ylim  (default 98 => pct=99)
%
% Figures saved as  <pairTag>_avg_direction_smooth.png  (or category)

if nargin<1, resultsRoot = 'results'; end
if nargin<2, win        = 5;          end
if nargin<3, method     = 'movmean';  end
if nargin<4, pct        = 99;         end          % 1–99 percentile

% ----- colour & pretty-name maps ---------------------------------------
clr.LIP = [0.85 0.15 0.15];
clr.FEF = [0.15 0.30 0.85];
clr.SC  = [0.15 0.60 0.15];
pretty  = @(tag) regexprep(tag, {'MLIP','MFEF','MSC'}, {'LIP','FEF','SC'});

pairDirs = dir(resultsRoot);
pairDirs = pairDirs([pairDirs.isdir] & ~startsWith({pairDirs.name},'.'));

% ---------- PASS 1 : global y-limits on **smoothed** data --------------
lim.direction = [ inf, -inf];
lim.category  = [ inf, -inf];

for p = 1:numel(pairDirs)
    pairTag = pairDirs(p).name;
    mats    = dir(fullfile(resultsRoot,pairTag,'session_mat','sess_*.mat'));
    if isempty(mats), continue; end

    for k = 1:numel(mats)
        s = load(fullfile(mats(k).folder,mats(k).name),'sess'); s=s.sess;
        if endsWith(pairTag,'direction')
            vals = smoothdata([s.dirAtoB(:); s.dirBtoA(:)], ...
                               method, win);
            hi   = prctile(vals, pct);
            lo   = prctile(vals, 100-pct);
            lim.direction(1) = min(lim.direction(1), lo);
            lim.direction(2) = max(lim.direction(2), hi);
        else
            vals = smoothdata([s.catAtoB(:); s.catBtoA(:)], ...
                               method, win);
            hi   = prctile(vals, pct);
            lo   = prctile(vals, 100-pct);
            lim.category(1)  = min(lim.category(1),  lo);
            lim.category(2)  = max(lim.category(2),  hi);
        end
    end
end

% ---------- PASS 2 : build figures -------------------------------------
for p = 1:numel(pairDirs)
    pairTag = pairDirs(p).name;
    toks    = split(pairTag,'_');
    srcTag  = toks{1};  dstTag = toks{2};  label = toks{3};
    mats    = dir(fullfile(resultsRoot,pairTag,'session_mat','sess_*.mat'));
    if isempty(mats), continue; end

    % gather & smooth each session trace --------------------------------
    for k = 1:numel(mats)
        s = load(fullfile(mats(k).folder,mats(k).name),'sess'); s=s.sess;
        if k==1, t=s.time_centers(:); A2B=[]; B2A=[]; end
        if strcmp(label,'direction')
            A2B=[A2B , smoothdata(s.dirAtoB(:) ,method,win)];
            B2A=[B2A , smoothdata(s.dirBtoA(:) ,method,win)];
        else
            A2B=[A2B , smoothdata(s.catAtoB(:) ,method,win)];
            B2A=[B2A , smoothdata(s.catBtoA(:) ,method,win)];
        end
    end

    mA = nanmean(A2B,2);  sA = nanstd(A2B,0,2);
    mB = nanmean(B2A,2);  sB = nanstd(B2A,0,2);
    ms = t*1e3;           N  = size(A2B,2);

    src = pretty(srcTag);  dst = pretty(dstTag);

    % plot --------------------------------------------------------------
    fig = figure('visible','off','Color','w'); hold on

    fill([ms; flipud(ms)], [mA+sA; flipud(mA-sA)], ...
         clr.(src),'FaceAlpha',0.18,'EdgeColor','none');
    plot(ms,mA,'-','Color',clr.(src),'LineWidth',2, ...
         'DisplayName',sprintf('%s \\rightarrow %s',src,dst));

    fill([ms; flipud(ms)], [mB+sB; flipud(mB-sB)], ...
         clr.(dst),'FaceAlpha',0.18,'EdgeColor','none');
    plot(ms,mB,'-','Color',clr.(dst),'LineWidth',2, ...
         'DisplayName',sprintf('%s \\rightarrow %s',dst,src));

    refline(0,0); grid on; box off
    xlabel('Time (ms)'); ylabel('FIT (bits)');
    title(sprintf('%s  (%s,  N=%d,  %s%d)', ...
          strrep(pairTag,'_','\_'), upper(label),N, method,win), ...
          'Interpreter','none');
    legend('Location','best');

    % unified y-scale
    if strcmp(label,'direction'), ylim(lim.direction);
    else                         ylim(lim.category);  end

    pngName = sprintf('avg_%s_smooth.png',label);
    saveas(fig, fullfile(resultsRoot,[pairTag '_' pngName]));
    close(fig);
    fprintf('[SMOOTH] %s  →  %s\n', pairTag, pngName);
end
end
