%rainf
% This script is designed to check if the plan is working

fid = fopen('basin_code.txt','r');
fom = '%d';
bid = fscanf(fid,fom);
fclose(fid);

tmp = load(['/shared/dondo/home/ml423/DA/Critical matrices/WORLD_events_v3.mat'],'WORLD_events');
WORLD_events = tmp.WORLD_events;

ixx = find(WORLD_events(bid,:)>0);
etot = length(ixx)-1;

fid = fopen('events_number.txt','w');
fprintf(fid,'%d',etot);
fclose(fid);

gid = WORLD_events(bid,1);
fid = fopen('basin_gauge_id.txt','w');
fprintf(fid,'%d',gid);
fclose(fid);

quit