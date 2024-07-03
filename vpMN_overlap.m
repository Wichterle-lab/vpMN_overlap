hix = 6;
mix = 3;
cix = 3; 

load ~/Dropbox/DATA/RNA/10X/Mouse-Human_final/01_compile/FINAL/Human.mat
load ~/Dropbox/DATA/RNA/10X/Mouse-Human_final/01_compile/FINAL/Mouse.mat
load ~/Dropbox/DATA/RNA/10X/Mouse-Human_final/01_compile/FINAL/Spinal.mat

load ~/Dropbox/DATA/RNA/10X/Mouse-Human_final/02_analysis/FINAL/human_umi.mat
load ~/Dropbox/DATA/RNA/10X/Mouse-Human_final/02_analysis/FINAL/mouse_umi.mat
load ~/Dropbox/DATA/RNA/10X/Mouse-Human_final/02_analysis/FINAL/human_genes.mat
load ~/Dropbox/DATA/RNA/10X/Mouse-Human_final/02_analysis/FINAL/mouse_genes.mat
load ~/Dropbox/DATA/RNA/10X/Mouse-Human_final/02_analysis/FINAL/human_cells.mat
load ~/Dropbox/DATA/RNA/10X/Mouse-Human_final/02_analysis/FINAL/mouse_cells.mat
htime = Human(hix).time;
mtime = Mouse(mix).time; 

% ho = unique(Human(hix).seurcid);
% mo = unique(Mouse(mix).seurcid);
% co = unique(Spinal(cix).seurcid);

ho = [6,8,2,5,3,1,4,0,9,7];
mo = [2,4,6,1,0,3,8,5,7,9];
co = [3,2,4,0,8,6,9,5,1,10,7,11];
load ~/Dropbox/MATLAB/FUNCS/viridi.mat

goi_h = {'LIN28A','SOX1','HES1','PAX6','OLIG1','OLIG2','NEUROG2','NEUROG1','LHX3','NKX2-2','SOX9','ISL1','MNX1','SLC18A3'};
goi_m = goi_h;

umap = Spinal(cix).umap; 
cca = Spinal(cix).seurcid;

ccl = cell(length(co),1);
for i = 1:length(ccl)
    ccl{i} = sprintf('C%d',co(i));
end

%%
h1 = Human(hix).seurcid;
h2 = Spinal(cix).seurcid;
c1 = Human(hix).idname;
c2 = Spinal(cix).idname;

h1 = h1(ismember(c1,c2)); 
c1 = c1(ismember(c1,c2));

h2 = h2(ismember(c2,c1)); 
c2 = c2(ismember(c2,c1));

hu1 = unique(h1);
hu2 = unique(cca);

if sum(ismember(c1,c2)) == length(c1)
    [c1,o] = sort(c1);
    h1 = h1(o);
    [c2,o] = sort(c2);
    h2 = h2(o);
    if isequal(c1,c2)
        [hm,y,x] = ccaclusteroverlap([h1,h2]); % row(y) = h1 values, col(x) = h2 values
    end
end
hm_h = hm;
save('heatmap_human.mat','hm_h');
% 

%%
m1 = Mouse(mix).seurcid;
m2 = Spinal(cix).seurcid;
c1 = Mouse(mix).idname;
c2 = Spinal(cix).idname;

m1 = m1(ismember(c1,c2)); 
c1 = c1(ismember(c1,c2));
m2 = m2(ismember(c2,c1)); 
c2 = c2(ismember(c2,c1));

if sum(ismember(c1,c2)) == length(c1)
    [c1,o] = sort(c1);
    m1 = m1(o);
    [c2,o] = sort(c2);
    m2 = m2(o);
    if isequal(c1,c2)
        [hm,y,x] = ccaclusteroverlap([m1,m2]); % row(y) = h1 values, col(x) = h2 values
    end
end

hm_m = hm;
save('heatmap_mouse.mat','hm_m');

%% which cca clusters are disproportionately represented by one species?
% temp = meta.Cellid;
% fh = find(contains(temp,'H'));
% fm = find(contains(temp,'M'));
% sp = cell(length(temp),1);hix = 6;


%% which cca clusters are disproportionately represented by one species?
% temp = meta.Cellid;
% fh = find(contains(temp,'H'));
% fm = find(contains(temp,'M'));
% sp = cell(length(temp),1);
% sp(fh) = {'human'};
% sp(fm) = {'mouse'};
% meta.Species = sp;
% 
% ucc = unique(cca);
% numhm = NaN(length(ucc), 2); 
% for i = 1:length(ucc)
%     f = cca == ucc(i);
%     temp = meta.Species;
%     numhm(i,1) = sum(ismember(temp(f),'human'));
%     numhm(i,2) = sum(ismember(temp(f),'mouse'));
% end

hm_h2 = hm_h./sum(hm_h,2);
hm_m2 = hm_m./sum(hm_m,2); 
phm = [max(hm_h2)', max(hm_m2)'];


%%
figure; 
subplot(2,1,1);

temp_h = hm_h;

hcl = cell(length(ho),1);
for i = 1:length(ho)
   hcl{i} = sprintf('H%d',ho(i)); 
end
heatmap(ccl,hcl,...
    temp_h((ho+1),(co+1)),'Colormap',viridi,...
    'XLabel','Common Clusters', 'YLabel', 'Human Clusters','GridVisible',true)

subplot(2,1,2);

temp_m = hm_m;
mcl = cell(length(mo),1);
for i = 1:length(mo)
   mcl{i} = sprintf('M%d',mo(i)); 
end
heatmap(ccl,mcl,...
    temp_m((mo+1),(co+1)),...
    'Colormap',viridi,...
    'XLabel','Common Clusters', 'YLabel', 'Mouse Clusters','GridVisible',true)

%%
figure; 
subplot(2,1,1);

temp_h = hm_h2;

hcl = cell(length(ho),1);
for i = 1:length(ho)
   hcl{i} = sprintf('H%d',ho(i)); 
end
heatmap(ccl,hcl,...
    temp_h((ho+1),(co+1)),'Colormap',viridi,'ColorLimits',[0.2,0.7],...
    'XLabel','Common Clusters', 'YLabel', 'Human Clusters','GridVisible',true)



subplot(2,1,2);
temp_m = hm_m2;
mcl = cell(length(mo),1);
for i = 1:length(mo)
   mcl{i} = sprintf('M%d',mo(i)); 
end
heatmap(ccl,mcl,...
    temp_m((mo+1),(co+1)),...
    'Colormap',viridi,'ColorLimits',[0.2,0.7],...
    'XLabel','Common Clusters', 'YLabel', 'Mouse Clusters','GridVisible',true)



%% fisher's exact test for common-species cluster relationships
h1 = Human(hix).seurcid;
h2 = Spinal(cix).seurcid;
c1 = Human(hix).idname;
c2 = Spinal(cix).idname;

h1 = h1(ismember(c1,c2)); 
c1 = c1(ismember(c1,c2));

h2 = h2(ismember(c2,c1)); 
c2 = c2(ismember(c2,c1));

if isequal(c1,c2)
    hu1 = unique(h1);
    hu2 = unique(cca);
    fet_h = NaN(length(hu1),length(hu2)); 
    for i = 1:length(hu1)
        f1 = h1 == hu1(i); 
        for j = 1:length(hu2)
            f2 = h2 == hu2(j); 
            x = NaN(2,2); 
            x(1,1) = sum((f1 == 1).*(f2 == 1));
            x(2,1) = sum((f1 == 0).*(f2 == 1));
            x(1,2) = sum((f1 == 1).*(f2 == 0));
            x(2,2) = sum((f1 == 0).*(f2 == 0));   
            [h,p,stats] = fishertest(x);
            fet_h(i,j) = -log(p)/log(10); 
        end
    end
end
figure; subplot(2,1,1)
heatmap(ccl,hcl,...
    fet_h((ho+1),(co+1)),...
    'Colormap',viridi,...
    'XLabel','Common Clusters', 'YLabel', 'Mouse Clusters','GridVisible',true)


m1 = Mouse(mix).seurcid;
m2 = Spinal(cix).seurcid;
c1 = Mouse(mix).idname;
c2 = Spinal(cix).idname;

m1 = m1(ismember(c1,c2)); 
c1 = c1(ismember(c1,c2));
m2 = m2(ismember(c2,c1)); 
c2 = c2(ismember(c2,c1));

if isequal(c1,c2)
    mu1 = unique(m1);
    mu2 = unique(cca);
    fet_m = NaN(length(mu1),length(mu2)); 
    for i = 1:length(mu1)
        f1 = m1 == mu1(i); 
        for j = 1:length(mu2)
            f2 = m2 == mu2(j); 
            x = NaN(2,2); 
            x(1,1) = sum((f1 == 1).*(f2 == 1));
            x(2,1) = sum((f1 == 0).*(f2 == 1));
            x(1,2) = sum((f1 == 1).*(f2 == 0));
            x(2,2) = sum((f1 == 0).*(f2 == 0));   
            [h,p,stats] = fishertest(x);
            fet_m(i,j) = -log(p)/log(10); 
        end
    end
end
subplot(2,1,2)
heatmap(ccl,mcl,...
    fet_m((mo+1),(co+1)),...
    'Colormap',viridi,...
    'XLabel','Common Clusters', 'YLabel', 'Mouse Clusters','GridVisible',true)


%% chi-square 
csq_h = NaN(size(hm_h,1), size(hm_m,1));
csq_m = NaN(size(hm_h,1), size(hm_m,1));
for i = 1:size(hm_h,1)
    v1 = hm_h(i,:);
    for j = 1:size(hm_m,1)
        v2 = hm_m(j,:);
        vsum  = v1 + v2; 
        vtot = sum(vsum); 
        v1e = vsum.*sum(v1)./vtot;
        v2e = vsum.*sum(v2)./vtot;
        csq = (v1 - v1e).^2./v1e; 
        csq = sum(csq(~isnan(csq)));
        csq_h(i,j) = csq; %chi2cdf(csq, 11); 
        csq = (v2 - v2e).^2./v2e; 
        csq = sum(csq(~isnan(csq)));
        csq_m(i,j) = csq;%chi2cdf(csq, 11); 
    end
end

temp = max(log(csq_h(:))) - log(csq_h);
figure; heatmap(mcl,hcl,temp((ho+1),(mo+1)),'Colormap',viridi,'ColorLimits',[1.5,3])

temp = log(csq_h);
e = size(viridi,1); v = e:-1:1;
figure; heatmap(mcl,hcl,temp((ho+1),(mo+1)),'Colormap',viridi(v,:),...
    'GridVisible',true,'ColorLimits', [3,4.8])

temp = log(csq_m); 
figure; heatmap(mcl,hcl,temp((ho+1),(mo+1)),'Colormap',viridi(v,:),...
    'GridVisible',true,'ColorLimits', [3,5])

temp = (log(csq_h)+log(csq_m))./2;
figure; heatmap(mcl,hcl,temp((ho+1),(mo+1)),'Colormap',viridi(v,:),...
    'GridVisible',true,'ColorLimits', [3,4.8])

temp = ((csq_h)+(csq_m))./2;
figure; heatmap(mcl,hcl,temp((ho+1),(mo+1)),'Colormap',viridi(v,:),...
    'GridVisible',true,'ColorLimits',[20,60])
