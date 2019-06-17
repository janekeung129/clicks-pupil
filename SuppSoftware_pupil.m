clear; clc;
load regs_pupil.mat

%% Figure 1
close all
clearvars -except regs

sigfun = @(b,x) b(1) + b(2)./ (1 + exp(-(x-b(3))/b(4)));
bin    = -14:2:14;
binabs = 0:2:14;
for sub = 1:length(regs)
    dn = sum(regs(sub).clicksLR);
    ch = regs(sub).choiceLR;
    rt = regs(sub).RT;
    for d = 1:length(bin)
        ind = dn==bin(d) & rt<3;
        pright(sub,d) = nanmean(ch(ind)==1);
        ntrial(sub,d) = sum(ind);
    end
    for d = 1:length(binabs)
        ind = abs(dn)==binabs(d) & rt<3;
        RT(sub,d)       = nanmean(rt(ind));
        ntrialrt(sub,d) = sum(ind);
    end
    prsub = pright(sub,~isnan(pright(sub,:)));
    rtsub = RT(sub,~isnan(RT(sub,:)));
    ch(ch==-1)=0;
    betapr(sub,:) = glmfit(dn,ch','binomial','link','logit');
    betart(sub,:) = glmfit(abs(dn(rt<3)),rt(rt<3));
    prightfit(sub,:) = glmval(betapr(sub,:)',[-14:.1:14]','logit');
    RTfit(sub,:) = glmval(betart(sub,:)',[0:2:14]','identity');
end

M = nanmean(pright);
S = nanstd(pright,[],1)/length(regs);
MM = nanmean(prightfit);
N = nanmean(ntrial);
M2 = nanmean(RT);
S2 = nanstd(RT,[],1)/length(regs);
MM2 = nanmean(RTfit);
N2 = nanmean(ntrialrt);

figure;
ax(1) = subplot(1,2,1); hold on;
p(1) = plot([-14:.1:14],MM);
sc(1) = scatter(bin,M,N/3);
errorbar(bin,M,S,'color','k','CapSize',1,'linestyle','none');
xlabel('clicks_{left}-clicks_{right}');
ylabel('p(choose left)');

ax(2) = subplot(1,2,2); hold on;
p(2) = plot(binabs,MM2);
sc(2) = scatter(binabs,M2,N2/3);
errorbar(binabs,M2,S2,'color','k','CapSize',1,'linestyle','none');
xlabel('|clicks_{left}-clicks_{right}|');
ylabel('reaction time (s)');

set(ax(1),'xlim',[-14.5 14.5],'ylim',[0 1],'ytick',[0 .5 1]);
set(ax(2),'xlim',[-.8 12.5],'ylim',[0 .8],'ytick',0:.2:1);
set(p,'color','k','linestyle','--');
set(sc,'MarkerEdgeColor','None','MarkerFaceColor',.65*ones(1,3));

%% Figure 2
close all
clearvars -except regs
nback = 5;
for sub = 1:length(regs)
    clicks = regs(sub).clicksLR';
    ch = regs(sub).choiceLR';
    dn = regs(sub).dN_all';
    ct = sign(dn);
    inds = repmat(1:length(ch)-nback,nback+1,1)+repmat([0:nback]',1,length(ch)-nback);
    ind_ct = inds(end,:)';
    ind_pt = inds(1:end-1,:)';
    y = ch(:);y(y==-1)= 0;
    x = [clicks(ind_ct,:) ch(ind_pt) ct(ind_pt)];
    y = y(ind_ct);
    betas(sub,:) = glmfit(x,y,'binomial','link','logit');
end
ind_dc = 2:21;
ind_ch = 22:22+nback-1;
ind_ct = 22+nback:22+nback+nback-1;

mat = {betas(:,ind_dc),betas(:,ind_ch),betas(:,ind_ct)};
dn = mean(mat{1},2);
ddn = mean(dn);
sdn = std(dn)/sqrt(size(mat{1},1));
M{1} = mean(mat{1}-repmat(dn,1,20));
S{1} = std(mat{1}-repmat(dn,1,20))/sqrt(size(mat{1},1));
for i = 2:3
    M{i} = mean(mat{i});
    S{i} = std(mat{i})/sqrt(size(mat{i},1));
end
side = mean(betas(:,1));
sside = std(betas(:,1))/sqrt(length(regs));

figure;
ax(1) = subplot(1,4,1); hold all;
errorbar(1,ddn,sdn,'k');
plot(ddn,'o','MarkerFaceColor','k','MarkerEdgeColor','none');
plot([0 21],[0 0],'k--');
ylabel('regression weight');

ax(4) = subplot(1,4,4); hold all;
plot(side,'o','MarkerFaceColor','k','MarkerEdgeColor','none');
errorbar(1,side,sside,'k');
plot([0 21],[0 0],'k--');

ax(2) = subplot(1,4,2); hold all;
plot([0 21],[0 0],'k--');
errorbar(1:length(M{1}),M{1},S{1},'k');
xlabel('click number');

ax(3) = subplot(1,4,3); hold all;
plot([0 6],[0 0],'k--');
errorbar(1:length(M{2}),M{2},S{2},'k','linestyle','none');
p(2) = plot(M{2},'k:','linewidth',2);
errorbar(1:length(M{3}),M{3},S{3},'k','linestyle','none');
p(3) = plot(M{3},'k--');
xlabel('previous trial number');

set(ax(1),'xlim',[0.5 1.5],'xtick',1,'xticklabel',{'SNR'});
set(ax(4),'xlim',[0.5 1.5],'xtick',1,'xticklabel',{'side bias'});
set(ax(2:4),'ytick',[]);
set(ax(2),'xlim',[0 21],'xtick', [1 5:5:20]);
set(ax(3),'xlim',[.5 5.5],'xtick', 1:5,'xticklabel',-5:-1);
legend(p([3 2]),{'RL','CK'},'location','southwest');legend boxoff
set(ax,'ylim',[-.35 .5]);

%% Figure 3
close all;
clearvars -except regs betas ind*
dn = mean(betas(:,2:21),2);
ch = betas(:,26);
ct = betas(:,31);
sc = {[ct,dn],[ch,dn],[ct,ch]};
ttl = {'SNR and RL','SNR and CK','RL and CK'};
xl  = {'RL','CK','RL'};
yl  = {'SNR','SNR','CK'};
for i = 1:length(sc)
    ax(i) = subplot(1,length(sc),i); hold all;
    scatter(sc{i}(:,1),sc{i}(:,2));
    lsline;
    title(ttl{i});
    ylabel([yl{i} ' regression weights']);
    xlabel([xl{i} ' regression weights']);
end
set(ax(1:2),'ylim',[min(dn(:))-.04 max(dn(:))+.04]);
set(ax(1),'xlim',[min(ct(:))-.04 max(ct(:))+.04]);
set(ax(2),'xlim',[min(ch(:))-.04 max(ch(:))+.04]);


%% Figure 4
close all
clearvars -except regs betas ind*

thr_trial = 60;
for sub = 1:length(regs)
    e = shiftdim(regs(sub).clicks_on_pd_z.D,2);
    e = e(:,31:120);
    keeptrial = sum(isnan(e),2)<thr_trial;
    erp(sub,:) = nanmean(e(keeptrial,:));
    ech_s{sub} = max(e(:,31:61),[],2) - min(e(:,31:61),[],2);
    ech(sub) = max(erp(sub,31:61)) - min(erp(sub,31:61));
end
E = ech;
ind1 = E>=prctile(E,50);
ind2 = E<prctile(E,50);

cblue = 'b';
cred = 'r';
colors = {'b','r'};

figure;
ax(4) = subplot(1,2,1); hold all;
merp  = nanmean(erp);
serp  = nanstd(erp)/sqrt(size(erp,1));
ph(1) = patch([1:length(merp) length(merp):-1:1],[merp+serp fliplr(merp-serp)],'k','EdgeColor','none');
p(7)  = plot(merp,'k');
plot([30 30],[-1 1],'k--');
plot([60 60],[-1 1],'k--');
xlabel('time (s)');
ylabel('pupil diameter (z-scored)');

ax(5) = subplot(1,2,2); hold all;
merp1 = nanmean(erp(ind1,:));
serp1 = nanstd(erp(ind1,:))/sqrt(sum(ind1));
merp2 = nanmean(erp(ind2,:));
serp2 = nanstd(erp(ind2,:))/sqrt(sum(ind2));
ph(3) = patch([1:length(merp2) length(merp2):-1:1],[merp2+serp2 fliplr(merp2-serp2)],[.9 .9 .9],'EdgeColor','none');
p(8)  = plot(merp2,'color',cred);
ph(2) = patch([1:length(merp1) length(merp1):-1:1],[merp1+serp1 fliplr(merp1-serp1)],[.9 .9 .9],'EdgeColor','none');
p(9)  = plot(merp1,'color',cblue);
plot([30 30],[-1 1],'k--');plot([60 60],[-1 1],'k--');
xlabel('time (s)');
legend(p([9 8]),{'high pupil change','low pupil change'},'location','northeast'); legend boxoff
set(ax(5),'ytick',[]);
set(ph(2),'facecolor',cblue);
set(ph(3),'facecolor',cred);
set(ph,'facealpha',.3);

figure;
mat_1 = {betas(ind1,ind_dc),betas(ind2,ind_dc)};
mat_2 = {betas(ind1,ind_ch),betas(ind1,ind_ct),betas(ind2,ind_ch),betas(ind2,ind_ct)};
for i = 1:2
    dn = mean(mat_1{i},2);
    ddn{i} = mean(dn);
    sdn{i} = std(dn)/sqrt(size(mat_1{1},1));
    M{i} = mean(mat_1{i}-repmat(dn,1,20));
    S{i} = std(mat_1{i}-repmat(dn,1,20))/sqrt(size(mat_1{i},1));
end
side{1} = mean(betas(ind1,1));
sside{1} = std(betas(ind1,1))/sqrt(sum(ind1));
side{2} = mean(betas(ind2,1));
sside{2} = std(betas(ind2,1))/sqrt(sum(ind2));
for i = 1:4
    Mp{i} = mean(mat_2{i});
    Sp{i} = std(mat_2{i})/sqrt(size(mat_2{i},1));
end
for i = 2:-1:1
    ax(1) = subplot(1,4,1); hold all;
    errorbar(1,[ddn{i}],[sdn{i}],'color',colors{i});
    plot([ddn{i}],'o','MarkerFaceColor',colors{i},'MarkerEdgeColor','none');
    plot([0 21],[0 0],'k--');
    ylabel('regression weight')

    ax(8) = subplot(1,4,4); hold all;
    errorbar(1,[side{i}],[sside{i}],'color',colors{i});
    plot([side{i}],'o','MarkerFaceColor',colors{i},'MarkerEdgeColor','none');
    plot([0 21],[0 0],'k--');
    
    ax(2) = subplot(1,4,2); hold all;
    plot([0 21],[0 0],'k--');
    errorbar(1:length(M{i}),M{i},S{i},'color',colors{i});
    p(i*3-2) = plot(M{i},'color',colors{i});
    xlabel('click position');
    
    ax(3) = subplot(1,4,3); hold all;
    plot([0 6],[0 0],'k--');
    errorbar(1:length(Mp{i*2-1}),Mp{i*2-1},Sp{i*2-1},'color',colors{i},'linestyle','none');
    p(i*3-1) = plot(Mp{i*2-1},':','color',colors{i});
    errorbar(1:length(Mp{i*2}),Mp{i*2},Sp{i*2},'color',colors{i},'linestyle','none');
    p(i*3) = plot(Mp{i*2},'--','color',colors{i});
    xlabel('previous trial position'); 
end
set(ax([1:3 8]),'ylim',[-.35 .5]);
set(ax([1 8]),'xlim',[0.5 1.5],'xtick',1);
set(ax(1),'xticklabel',{'snr'});
set(ax(8),'xticklabel',{'side bias'});
set(ax([2 3 8]),'ytick',[]);
set(ax(2),'xlim',[0 21],'xtick', [1 5:5:20]);
set(ax(3),'xlim',[.5 5.5],'xtick', 1:5,'xticklabel',-5:-1);
set(ax(4:5),'xlim',[.5 90.5]);
set(ax(4:5),'ylim',[-.27 .37],'ytick',-1:.2:1,...
    'xtick', [1 30:30:90],'xticklabel',-1:2);

%% Figure 5
close all
clearvars -except regs e*

for sub = 1:length(regs)
    clicks = regs(sub).clicksLR';
    ch = regs(sub).choiceLR';
    dn = regs(sub).dN_all';
    ct = sign(dn);
    e  = ech_s{sub};
    len = 1:min(size(clicks,1),length(e));
    nback = 1;
    inds = repmat(1:length(len)-nback,nback+1,1)+repmat([0:nback]',1,length(len)-nback);
    ind_ct = inds(end,:)';
    ind_pt = inds(1:end-1,:)';
    x = [clicks(ind_ct,:) ch(ind_pt) ct(ind_pt) dn(ind_ct).*e(ind_ct) ch(ind_pt).*e(ind_ct) ct(ind_pt).*e(ind_ct) e(ind_ct)];
    y = ch(ind_ct);y(y==-1)= 0;   
    betast(sub,:) = glmfit(x,y,'binomial','link','logit');
end

figure;
hold all;
M = mean(betast(:,[27 24 26 25]));
S = std(betast(:,[27 24 26 25]))./sqrt(length(regs));
errorbar(1:4,M,S,'color','k','linestyle','none');
plot(M,'ok','MarkerFaceColor','k');
plot([0 5],[0 0],'k--');
ylabel('regression weight')
set(gca,'xlim',[.5 4.5],'xtick',1:4,'ylim',[-0.05 0.05],'xticklabel',...
    {'pupil','\Delta click x pupil','previous correct x pupil','previous choice x pupil'});

%% Figure 6
close all
clearvars -except regs e*
load fits.mat

for s = 1:length(list)
    sub = list(s);
    X = regs(sub).clicksLR';
    y = regs(sub).choiceLR';y(y==-1)= 0;
    LLs = LLssub{s};
    pr = zeros(length(y),1);
    pr(y==0) = 1-exp(LLs(y==0));
    pr(y==1) = exp(LLs(y==1));
    pl = 1-pr;
    chddm{s} = double(~binornd(1,pl));
    betasddm(s,:) = glmfit(X,chddm{s},'binomial','link','logit');
    betas(s,:) = glmfit(X,y,'binomial','link','logit');
end
ind_dc = 2:21;

figure;
colors = {'b','r'}; ls = {'-','--'};
M = {mean(betas(:,ind_dc),1),mean(betasddm(:,ind_dc),1)};
S = {std(betas(:,ind_dc),1)/sqrt(length(list)),std(betasddm(:,ind_dc),1)};
ax(1) = subplot(1,3,1); hold on;
for i = 1:2
    e(i) = patch([1:20 20:-1:1],[S{i}+M{i} fliplr(M{i}-S{i})],colors{i});
    p(i) = plot(M{i},ls{i},'color',colors{i});
end
legend(p,'human','DDM','Location','southeast'); legend boxoff
set(e,'facealpha',0.3,'edgecolor','none');
ylabel('regression weight');
xlabel('click position');
set(gca,'xlim',[0 21],'xtick', [1 5:5:20],'ylim',[.0 .61],'ytick',0:.2:.8);

ax(2) = subplot(1,3,2); hold on;
nsub = length(list);
bin = -14:2:14;
for s = 1:length(list)
    sub = list(s);
    clear ch
    dn = sum(regs(sub).clicksLR);
    ch(1,:) = regs(sub).choiceLR;
    ch(2,:) = chddm{s}';
    rt = regs(sub).RT;
    for d = 1:length(bin)
        ind = dn==bin(d) & rt<3;
        for t = 1:2, pright(s,d,t) = nanmean(ch(t,ind)==1);end;
        ntrial(s,d) = sum(ind);
    end
end
for t = 1:2
    pr = nanmean(squeeze(pright(:,:,t)));
    sd = nanstd(squeeze(pright(:,:,t)))./sqrt(length(regs));
    sc(t) = scatter(1:15,pr,nanmean(ntrial)./3,'markerfacecolor','none','markeredgecolor',colors{t});
    errbar(1:15,pr,sd,'color',colors{t});
end
set(ax(2),'xlim',[.5 15.5],'xtick',[3 8 13],'xticklabel',{'-10','0','10'},'ylim',[-.02 1.02],'ytick',[0 .5 1]);
xlabel('clicks_{left}-clicks_{right}');
ylabel('p(choose left)');
legend(sc,'human','DDM','Location','southeast');legend boxoff

ax(3) = subplot(1,3,3); hold on;
sc(2) = scatter(paramddm(:,5),ech',15,'markerfacecolor',colors{1},'markeredgecolor','none');
lsline;
ylabel('pupil change (z-scored)');
xlabel('bound');
ylim = get(ax(3),'YLim');
set(ax(3),'xlim',[5 20],'YLim',[ylim(1)-diff(ylim)*.005 ylim(2)+diff(ylim)*.005]);
