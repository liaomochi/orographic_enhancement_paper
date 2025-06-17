%rainf
% This script is designed to check if the plan is working
fid = fopen('current_basin.txt','r');
fom = '%d';
basinid = fscanf(fid,fom);
fclose(fid); 

tmp = load(['/shared/dondo/home/ml423/DA/Critical matrices/world_selected_gauge_events_allinfo_v3.mat'],'gauge_events_first_last');
Basins_sum = tmp.gauge_events_first_last;clear tmp;
tmp = load(['/shared/dondo/home/ml423/DA/Critical matrices/WORLD_events_v3.mat'],'WORLD_events');
WORLD_events = tmp.WORLD_events;

gid = WORLD_events(basinid,1);

gloc = find(Basins_sum(:,1)==gid);

fid = fopen('datnm.txt','r');
fom = '%d';
event = fscanf(fid,fom);
fclose(fid);

nmm = num2str(event);
tww = 3;


tmp = load(['/shared/dondo/home/ml423/world_mts/obs_events/events/obs_' num2str(gid) '_event' num2str(event) '.mat']);
tmp = tmp.obsnow;
obsnow = tmp;clear tmp

obsnan=find(isnan(obsnow));
if isempty(obsnan)
else
    for na = 1:length(obsnan)
        obsnow(obsnan(na)) = 0.5*(obsnow(obsnan(na)-1)+obsnow(obsnan(na)+1));
    end
end

obs2 = interp1(1:31,obsnow,[1:0.33333:31]);
obsnow = obs2(1:90);


ffnm=['/shared/dondo/home/ml423/world_mts_outputs/Basin' num2str(gid) '_output_900m1hr/' nmm '_SW/' nmm '/streamflowarea'];
fid=fopen(ffnm,'rb','ieee-le');
tmpdata=fread(fid,inf,'single');
NODAc2o1=reshape(tmpdata,Basins_sum(gloc,7),Basins_sum(gloc,6),720);
fclose(fid);

for k = 1:720
    flownoda3(k) = NODAc2o1(Basins_sum(gloc,9),Basins_sum(gloc,8),k);
end

no3 = flownoda3;
simu = no3(1:8:720);


rr = corrcoef(simu',obsnow);
rstar = rr(1,2);
obsstd = std(obsnow);
simustd = std(simu);
obsmean = mean(obsnow);
simumean = mean(simu);
KGE = 1-sqrt((rstar-1)^2+(simustd/obsstd-1)^2+(simumean/obsmean-1)^2);


if KGE<=0.95
    moreitr = 1;
else
    moreitr = 0;
end

fid = fopen('moreitr.txt','w');
fprintf(fid,'%d',moreitr);
fclose(fid);

quit