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

fid = fopen('utchh.txt','r');
fom = '%d';
utch = fscanf(fid,fom);
fclose(fid);

tmp = load(['/shared/dondo/home/ml423/DA/Critical matrices/World_Basins/distance/Basin' num2str(gid) 'dis_stream.mat']);
tmp = tmp.dis_stream;
dis_stream = tmp;
tmp = load(['/shared/dondo/home/ml423/DA/Critical matrices/World_Basins/distance/Basin' num2str(gid) 'dis_outlet.mat']);
tmp = tmp.dis_outlet;
dis_outlet = tmp;


tmp = load(['/shared/dondo/home/ml423/DA/Critical matrices/World_Basins/inttind/Basin' num2str(gid) '_ind.mat'],'ind');
ind = tmp.ind;
tmp = load(['/shared/dondo/home/ml423/DA/Critical matrices/World_Basins/inttind/Basin' num2str(gid) '_intt.mat'],'intt');
intt = tmp.intt;



ffnm=['/shared/dondo/home/ml423/world_mts/Basin' num2str(gid) '_fixed_900m/soilparameters.ascii'];
soi = dlmread(ffnm);

wcs1 = reshape(soi(:,5),Basins_sum(gloc,7),Basins_sum(gloc,6));
wcs2 = reshape(soi(:,6),Basins_sum(gloc,7),Basins_sum(gloc,6));
wcs3 = reshape(soi(:,7),Basins_sum(gloc,7),Basins_sum(gloc,6));
wcs4 = reshape(soi(:,8),Basins_sum(gloc,7),Basins_sum(gloc,6));

wit1 = reshape(soi(:,13),Basins_sum(gloc,7),Basins_sum(gloc,6));
wit2 = reshape(soi(:,14),Basins_sum(gloc,7),Basins_sum(gloc,6));
wit3 = reshape(soi(:,15),Basins_sum(gloc,7),Basins_sum(gloc,6));
wit4 = reshape(soi(:,16),Basins_sum(gloc,7),Basins_sum(gloc,6));

ffnm=['/shared/dondo/home/ml423/world_mts/Basin' num2str(gid) '_fixed_900m/soildepth.ascii'];
soi = dlmread(ffnm);
dsm = reshape(soi(:,1),Basins_sum(gloc,7),Basins_sum(gloc,6));
dp = reshape(soi(:,2),Basins_sum(gloc,7),Basins_sum(gloc,6));
ds0 = 0.1*ones(Basins_sum(gloc,7),Basins_sum(gloc,6));

wsmax = wcs1.*ds0;
wsmin = wit1.*ds0;
wsmmax = wcs2.*dsm;
wsmmin = wit2.*dsm;
wsdmax = wcs3.*dp;
wsdmin = wit3.*dp;


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
fclose(fid);

fid = fopen('flow_dis_str.txt','r');
fom = '%d';
critnum = fscanf(fid,fom);
fclose(fid);

fid = fopen('pre_rise_pt.txt','r');
fom = '%d';
rains_bd = fscanf(fid,fom);
fclose(fid);


% load input initial conditions
ffnm=['/shared/dondo/home/ml423/world_mts_outputs/Basin' num2str(gid) '_output_900m1hr_IC/' nmm 'c/TrueIC/laststep_ws'];
fid=fopen(ffnm,'rb','ieee-le');
tmpdata=fread(fid,inf,'single');
ws=reshape(tmpdata,Basins_sum(gloc,7),Basins_sum(gloc,6));
ffnm=['/shared/dondo/home/ml423/world_mts_outputs/Basin' num2str(gid) '_output_900m1hr_IC/' nmm 'c/TrueIC/laststep_wsm'];
fid=fopen(ffnm,'rb','ieee-le');
tmpdata=fread(fid,inf,'single');
wsm=reshape(tmpdata,Basins_sum(gloc,7),Basins_sum(gloc,6));
ffnm=['/shared/dondo/home/ml423/world_mts_outputs/Basin' num2str(gid) '_output_900m1hr_IC/' nmm 'c/TrueIC/laststep_wsd'];
fid=fopen(ffnm,'rb','ieee-le');
tmpdata=fread(fid,inf,'single');
wsd=reshape(tmpdata,Basins_sum(gloc,7),Basins_sum(gloc,6));


wso = ws;
wsmo = wsm;
% done input initial conditions

wron = load(['/shared/dondo/home/ml423/world_mts_runs/IRC1_sup/wrontest5MMORIw' num2str(ceil(p/tww)) num2str(mod(filnm,10)) 'cc' num2str(filnm) '_' num2str(gid) '_' num2str(utch,'%2.2d') '.mat']);
wron = wron.wron;

% pbb2 = load(['/shared/dondo/home/ml423/world_mts_runs/IRC1_sup/pbb2test5MMORIw' num2str(ceil(p/tww)) num2str(mod(filnm,10)) 'cc' num2str(filnm) '_' num2str(gid) '_' num2str(utch,'%2.2d') '.mat'],'pbb2');
% pbb2 = pbb2.pbb2;


workingwindow = 3;
modd = 2;

fid = fopen('inflection_pt.txt','r');% here is very different###########!!!!!!!!!
fom = '%d';
ck = fscanf(fid,fom);
fclose(fid);
% now the ck is the end point of the window 4, the early recession window,
% hydrograph after the early recession window, will be controlled by the SM
% UPSTREAM to the tracked pixels.

disp(['=============== inflection is point ' num2str(ck)])


pao=1;


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

fid = fopen('loop2.txt','r');
fom = '%d';
itr = fscanf(fid,fom);
fclose(fid);

% fid = fopen('loop3.txt','r');
% fom = '%d';
% lp3 = fscanf(fid,fom);
% fclose(fid);


% mode 1 all 10 replicates
% obs2mat = load('obs2mat.mat');
% obs2mat = obs2mat.obs2mat;
% mode 2 coarse hydrograph
% obsdata = load('b2agg10.mat');
% obsdata = obsdata.obsdatav4;
% original observation


% obsnow = obs2mat(:,lp3);
% disp(['$$$$$$$$$$$ loop3 value is ' num2str(lp3)])


ak = 4+(ck-1)*8;

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

percx = (simu(ck)-obsnow(ck))./obsnow(ck);
perc = abs(simu(ck)-obsnow(ck))./obsnow(ck);

dif11 = simu(ck-1)-obsnow(ck-1);
dif22 = simu(ck)-obsnow(ck);

difv = (dif11+dif22)*28800/2;% 8 hour is 28800 seconds
disp(['$$$ water diff after sm change is ' num2str(difv)])

tmp = 0;
tmp1 = 0;

tmp2 = mean(obsnow);
for i = 1:ck
    tmp = tmp + (simu(i)-obsnow(i))^2;
    tmp1 = tmp1 + (tmp2-obsnow(i))^2;
end
NSE = 1-tmp./tmp1;




cdy = 2;

if modd == 2

    %dif1 = simu(90)-obsnow(90);
    %dif2 = simu(ck)-obsnow(ck);

    %rdif = (dif1+dif2)*(90-ck)*15*60;

    veryslow = min(ck+6,90);

    rdif = (sum(simu(ck:veryslow))-sum(obsnow(ck:veryslow)))*8*3600;


    disp(['$$run point is ' num2str(ck)])
    disp(['$$ simu value is ' num2str(simu(ck))])
    
    disp(['$$ water diff value is ' num2str(rdif)])
    %disp(['$$$$$$$$$$$ water diff value is ' num2str(rdif2)])
    %disp(['$$$$$$$$$$$ water diff value is ' num2str(rdif3)])

    %end
    % mod 1
    %rainshare = rdif.*(pbb2(1:bdy)./sum(pbb2(1:bdy)));
    % mod 2
    %rainshare = rdif.*ones(bdy,1)./bdy;
    % one can consider to add iterative processes

    disp(['=============== NSE is ' num2str(NSE)])




    %rainshare = rdif(part).*ones((bdy-cdy+1),1)./(bdy-cdy+1);
    %rainshare = rdif(part).*(pbb2(cdy:bdy)./sum(pbb2(cdy:bdy)));
    rainshare = rdif;%.*[1/5*ones(1,48-cdy+1)./(48-cdy+1),4/5*(pbb2(49:bdy)./sum(pbb2(49:bdy)))];


    for jh = cdy%ak:-1:cdy
        %
        % disp(['$$$$$ state value is ' num2str(ak-(jh-1))])
        %gro = sum(sum(wron{ak,ak-(jh-1)}(:,:)));
        %          if (ak-(jh-1))>48%54%54 very important
        wp = wron{ak,ak-(jh-1)};
        if rains_bd==1
            wp = wron{ak,rains_bd+1}+wron{ak,rains_bd};%3670
        elseif rains_bd==0
            wp = wron{ak,rains_bd+1};
        else
            wp = wron{ak,rains_bd+1}+wron{ak,rains_bd-1}+wron{ak,rains_bd};%3670
        end

        % wp = zeros(Basins_sum(gloc,7),Basins_sum(gloc,6));
        % 
        % for interpo = -4:4
        %     for witr = -4:4
        %         if (rains_bd+witr)<ak
        %             wp = wp + wron{ak+interpo,rains_bd+witr+interpo};
        %             disp(['(ak, ) pairs are' num2str(ak,rains_bd+witr)])
        % 
        %         else
        %             disp(['temporarily out of boundary (ak, ) pairs are' num2str(ak) ' ' num2str(rains_bd+witr)])
        %             % slightly adjust ak to be a later timing than rains_bd
        %         end
        %     end
        % end

        [rowi,colj] = find(wp>=1);
        if ~isempty(rowi)
            % This is to find the earliest launched particles,
            % backtracking to the earliest will make sure the maximum
            % area of the impact area, that IC can make an impact on
            % the hydrograph.
            eligi_pt = find(wp>=1);
            dis_tan = dis_outlet(eligi_pt);
            dis_str = dis_stream(eligi_pt);
            if ismember(0, dis_str)&&length(eligi_pt)>=5
                distance_thre = min(dis_outlet(eligi_pt(find(dis_str==0))));
                % this threshold is the first occurence of stream pixel
                % that reached to the hydrograph curvature changing
                % point, so any pixel with distance longer than this
                % one needs to be SM modified during slow recession.
                break
            else
                distance_thre = 0.5*max(dis_outlet(:));% added 08=03-2024 to avoid the no definition of distance_thre
            end
        else
            distance_thre = 0.5*max(dis_outlet(:));% added 08=03-2024 to avoid the no definition of distance_thre
        end
    end

    %{
        if isempty(rowi)
            disp('SM adjustment is Jumped')
            
        else
            allpts = [];
            neigh8 = [i-1,j;i+1,j;i-1,j-1;i,j-1;i+1,j-1;i-1,j+1;i,j+1;i+1,j+1];
            ups = [rowi,colj];
            for itrr = 1:1000
                vps = ups;
                uss = size(vps,1);
                ups = [];
                allpts = [allpts;vps];

                for s = 1:uss
                    i = vps(s,1);
                    j = vps(s,2);

                    if fdr(i-1,j)==2
                        ups = [ups;[i-1,j]];
                    end
                    if fdr(i+1,j)==32
                        ups = [ups;[i+1,j]];
                    end
                    if fdr(i-1,j-1)==4
                        ups = [ups;[i-1,j-1]];
                    end
                    if fdr(i,j-1)==8
                        ups = [ups;[i,j-1]];
                    end
                    if fdr(i+1,j-1)==16
                        ups = [ups;[i+1,j-1]];
                    end
                    if fdr(i-1,j+1)==1
                        ups = [ups;[i-1,j+1]];
                    end
                    if fdr(i,j+1)==128
                        ups = [ups;[i,j+1]];
                    end
                    if fdr(i+1,j+1)==64
                        ups = [ups;[i+1,j+1]];
                    end
                end
                if isempty(ups)
                    break
                end

            end


            unimpactarea = zeros(76,84);
            for i = 1:size(allpts,1)
                unimpactarea(allpts(i,1),allpts(i,2))=1;
            end



            final_corr_area = zeros(76,84);
            final_corr_area(find(unimpactarea>=1&dis_stream==0)) = 1;
    %}


    % above is to find out impacted area based on flow time and
    % dis_stream
    final_corr_area = zeros(Basins_sum(gloc,7),Basins_sum(gloc,6));
    final_corr_area(find(dis_stream==critnum&dis_outlet>=distance_thre)) = 1;

    tmpp1 = find(final_corr_area>0);
    [l1,l2] = find(final_corr_area>0);
    %tmpp2 = find(rainp(tmpp1)>0);
    %gro = length(tmpp2);

    tmpp2 = 1:length(tmpp1);
    gro = length(tmpp1);

    sumsm = sum(wsmo(tmpp1(tmpp2))-wsmmin(tmpp1(tmpp2)));% careful using second layer!
    disp(['$$ slow SM correction is ' num2str(rainshare)])
    disp(['$$ impacted number pixel is ' num2str(length(tmpp1))])

    if sumsm>0
        for ju = 1:gro
            correction = rainshare/900/900*(wsmo(tmpp1(tmpp2(ju)))-wsmmin(tmpp1(tmpp2(ju))))./sumsm;
            %correction = rainshare(ak-(jh-1))/250/250*1000*rainp(tmpp1(tmpp2(ju)))./sumrain;
            wsm(l1(tmpp2(ju)),l2(tmpp2(ju))) = wsm(l1(tmpp2(ju)),l2(tmpp2(ju)))-correction;
            ws(l1(tmpp2(ju)),l2(tmpp2(ju))) = ws(l1(tmpp2(ju)),l2(tmpp2(ju)))-correction;
            wsd(l1(tmpp2(ju)),l2(tmpp2(ju))) = wsd(l1(tmpp2(ju)),l2(tmpp2(ju)))-correction;

        end

    end

    wsm(find(wsm<wsmmin)) = wsmmin(find(wsm<wsmmin));
    ws(find(ws<wsmin)) = wsmin(find(ws<wsmin));
    wsd(find(wsd<wsdmin)) = wsdmin(find(wsd<wsdmin));

    wsm(find(wsm>wsmmax)) = wsmmax(find(wsm>wsmmax));
    ws(find(ws>wsmax)) = wsmax(find(ws>wsmax));
    wsd(find(wsd>wsdmax)) = wsdmax(find(wsd>wsdmax));
end



%end

ffnm=['/shared/dondo/home/ml423/world_mts_outputs/Basin' num2str(gid) '_output_900m1hr_IC/' nmm 'c/TrueIC/laststep_ws'];
fid = fopen(ffnm,'wb');
fwrite(fid,ws,'single');
fclose(fid);
ffnm=['/shared/dondo/home/ml423/world_mts_outputs/Basin' num2str(gid) '_output_900m1hr_IC/' nmm 'c/TrueIC/laststep_wsm'];
fid = fopen(ffnm,'wb');
fwrite(fid,wsm,'single');
fclose(fid);
ffnm=['/shared/dondo/home/ml423/world_mts_outputs/Basin' num2str(gid) '_output_900m1hr_IC/' nmm 'c/TrueIC/laststep_wsd'];
fid = fopen(ffnm,'wb');
fwrite(fid,wsd,'single');
fclose(fid);

clear all;
%

quit
