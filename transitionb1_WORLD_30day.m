tmp = load(['/shared/dondo/home/ml423/DA/Critical matrices/world_selected_gauge_events_allinfo_v3.mat'],'gauge_events_first_last');
Basins_sum = tmp.gauge_events_first_last;clear tmp;

fid = fopen('current_basin.txt','r');
fom = '%d';
basinid = fscanf(fid,fom);
fclose(fid);

tmp = load(['/shared/dondo/home/ml423/DA/Critical matrices/WORLD_events_v3.mat'],'WORLD_events');
WORLD_events = tmp.WORLD_events;

gid = WORLD_events(basinid,1);

gloc = find(Basins_sum(:,1)==gid);

fid = fopen('utchh.txt','r');
fom = '%d';
utch = fscanf(fid,fom);
fclose(fid);

fid = fopen('datnm.txt','r');
fom = '%d';
event = fscanf(fid,fom);
fclose(fid);

nmm = num2str(event);
tww = 3;

fid = fopen('bl.txt','r');
fom = '%d';
p = fscanf(fid,fom);
fclose(fid);


fid = fopen('filenum.txt','r');
fom = '%d';
filnm = fscanf(fid,fom);
% ===== above is add for all

rainfallb1 = load('rainfallb1.mat');
rainfallb1 = rainfallb1.rainfallb1;

xtmpbb1 = rainfallb1;
save(['xtmpbb1.mat'],'xtmpbb1')

% corr = {};
% for i = 1:288
%    corr{i} = zeros(76,84); 
% end


%strr = ['corrtest5MMORIw' num2str(ceil(p/tww)) num2str(mod(filnm,10)) 'cc' nmm 'narrowwide'];
%save(['/shared/dondo/home/ml423/HMIOP/NoDAre/New era/' strr '.mat'],'corr')

%ounm = {'Precipitation_IMERGEnewb111'};
%ounm = {'Precipitation_St4dbknewb111'};
 ounm = {'Precipitation_ERA5_land_oriunit'};




for i = 1  %:46 

    datestr = nmm;

    out_dir = ['/shared/dondo/home/ml423/world_mts/Basin' num2str(gid) '_input_ERA5_900m_tmp/' datestr '/'];
    %out_dir = ['/shared/dondo/home/ml423/HMIOP/NoDAre/LS/Basin01_input_HistMODIS_250m5min/' datestr '/'];
    % out_dir = ['/shared/dondo/home/ml423/DA/LBK/LS/Basin01_input_WRF-AF_250m5min_Perturb/' datestr '/'];
    dataf = [];
    for j = 1:720
        datas = rainfallb1{j}; % can add stuff
        dataf(:,:,j) = datas;
    end
    dataf(dataf<0) = 0;
    fnm = [out_dir,ounm{1}];
    raindata = reshape(dataf,[Basins_sum(gloc,7),Basins_sum(gloc,6),720]);
    fid = fopen(fnm,'wb');
    fwrite(fid,raindata,'single');
    fclose(fid);
    fid = fopen([fnm,'_sta.txt'],'w');
    fprintf(fid,'Min    = %.4f\n',min(raindata(:)));
    fprintf(fid,'Max    = %.4f\n',max(raindata(:)));
    fprintf(fid,'Mean    = %.4f\n',mean(raindata(:)));
    fprintf(fid,'Std    = %.4f\n',std(raindata(:)));
    fclose(fid); 
end
quit