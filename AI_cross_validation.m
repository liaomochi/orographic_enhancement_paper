%% [AI Key!]Alps prediction using dis_outlet and dis_stream, use basin P as a mean
addpath('/shared/dondo/home/ml423/DA/slanCM/')
clear kgeori_median kge_median d_area real_event select_event
% exclude basins that are not in the Alps and failed the quality control (e.g. drainage area too big). 
for i = [1:129 132:133 466:480 482:504 507:522]
    eqcloc = find(~isnan(eqcontrol(i,:)));
    kgeori_median(i) = nanmedian(kge_orim(i,eqcloc)); 
    kge_median(i) = nanmedian(kge_optm(i,eqcloc)); 
    d_area(i,1) = gauge_events_first_last(find(gauge_events_first_last(:,1)==WORLD_events(i,1)),5);
    real_event(i) = length(find(kge_optm(i,eqcloc)>-999));
    select_event(i) = length(find(WORLD_events(i,eqcloc)>-999));
end

figure
bar([kgeori_median',kge_median'])
figure
bar(kgeori_median)
figure
bar(kge_median)

% exclude basins from different regions (outside the Alps) in the Europe (e.g.Pyrenneys) below!
problem_basins = [1:23 94 466 467 468 470 471 495 496 497 498 499 500];

figure
bar(d_area)
figure
scatter(select_event,real_event,30,d_area,'filled')


tmp = load(['/shared/dondo/home/ml423/DA/Critical matrices/world_selected_gauge_events_allinfo_v3.mat'],'gauge_events_first_last');
Basins_sum = tmp.gauge_events_first_last;clear tmp;
tmp = load(['/shared/dondo/home/ml423/DA/Critical matrices/WORLD_events_v3.mat'],'WORLD_events');
WORLD_events = tmp.WORLD_events;clear tmp;

% flowtype = {'streamflowarea','interflowarea','overlandflowarea','baseflowarea'};
% ftp = 1;

d_s_sum = {};
d_o_sum = {};
for bid = [1:129 132:133 466:480 482:504 507:522]

    gid = WORLD_events(bid,1);
    
    gloc = find(Basins_sum(:,1)==gid);
    
    r = Basins_sum(gloc,6);
    c = Basins_sum(gloc,7);
    rg = Basins_sum(gloc,8);
    cg = Basins_sum(gloc,9);
    area = Basins_sum(gloc,5);
    
    tmp = load(['/shared/dondo/home/ml423/DA/Critical matrices/World_Basins/inttind/Basin' num2str(gid) '_ind.mat'],'ind');
    ind = tmp.ind;
    tmp = load(['/shared/dondo/home/ml423/DA/Critical matrices/World_Basins/inttind/Basin' num2str(gid) '_intt.mat'],'intt');
    intt = tmp.intt;       

    d_s = load(['/shared/dondo/home/ml423/DA/Critical matrices/World_Basins/distance/Basin' num2str(gid) 'dis_stream.mat']);
    d_s = d_s.dis_stream;
    d_o = load(['/shared/dondo/home/ml423/DA/Critical matrices/World_Basins/distance/Basin' num2str(gid) 'dis_outlet.mat']);
    d_o = d_o.dis_outlet;
    d_s(ind) = nan;d_o(ind) = nan;
    d_s = d_s+1;%dis_stream is set to be 1 for minimum
    %d_s = d_s/max(d_s(:));
    d_o = d_o/max(d_o(:));
    
    d_s_sum{bid} = d_s;
    d_o_sum{bid} = d_o;
end

% look at dis_Stream
h_sum = {};
hsummax = [];
figure
for i = [1:129 132:133 466:480 482:504 507:522]
h1 = histogram(d_s_sum{i},'Normalization','probability');
h11 = h1.Values;
h_sum{i} = h11;
[~,maxloc] = max(h11);
hsummax(i) = maxloc;
hold on
end
outlie = find(hsummax==3)% the minority here that peaks at dis_stream=3

figure
for i = [1:129 132:133 466:480 482:504 507:522]
    plot(h_sum{i})
    hold on
end
title('dis-stream')

figure
for i = outlie
    plot(h_sum{i})
    hold on
end
% 6 outliers that dont peak at dis_stream=1
for i = outlie
    gid = WORLD_events(i,1);
    gloc = find(Basins_sum(:,1)==gid);
    
    r = Basins_sum(gloc,6);
    c = Basins_sum(gloc,7);
    tmp = load(['/shared/dondo/home/ml423/DA/Critical matrices/World_Basins/inttind/Basin' num2str(gid) '_ind.mat'],'ind');
    ind = tmp.ind;
    ffnm=['/shared/dondo/home/ml423/world_mts/Basin' num2str(gid) '_fixed_900m/dem.bin'];
    fid=fopen(ffnm,'rb','ieee-le');
    tmpdata=fread(fid,inf,'single');
    fclose(fid);
    dem = reshape(tmpdata,c,r);
    dem(ind) = nan;
    demf = dem';
    colo = slanCM('terrain',80);
    colo(1,:) = [0.8 0.8 0.8];
    figure
    imagesc(demf)
    colormap(colo)
end

% look at dis_outlet
h_oum = {};
figure
for i = 1:133
h1 = histogram(d_o_sum{i},'Normalization','probability');
h11 = h1.Values;
h12 = h1.BinEdges;
h_oum{i} = [h11;h12(1:end-1)];
 max(d_o_sum{i}(:))
hold on
end

figure
for i = 1:133
    plot(h_oum{i}(2,:),h_oum{i}(1,:))
    hold on
end

% calculate slopes
slopesum = {};
slopemeansum = [];
demsum = [];
for bid = [1:129 132:133 466:480 482:504 507:522]
    gid = WORLD_events(bid,1);
    gloc = find(Basins_sum(:,1)==gid);
    
    r = Basins_sum(gloc,6);
    c = Basins_sum(gloc,7);
    tmp = load(['/shared/dondo/home/ml423/DA/Critical matrices/World_Basins/inttind/Basin' num2str(gid) '_ind.mat'],'ind');
    ind = tmp.ind;
    tmp = load(['/shared/dondo/home/ml423/DA/Critical matrices/World_Basins/inttind/Basin' num2str(gid) '_intt.mat'],'intt');
    intt = tmp.intt; 
    ffnm=['/shared/dondo/home/ml423/world_mts/Basin' num2str(gid) '_fixed_900m/dem.bin'];
    fid=fopen(ffnm,'rb','ieee-le');
    tmpdata=fread(fid,inf,'single');
    fclose(fid);
    dem = reshape(tmpdata,c,r);
    %dem(ind) = nan;
    demsum(bid) = mean(dem(intt));

    ffnm=['/shared/dondo/home/ml423/world_mts/Basin' num2str(gid) '_fixed_900m/fdr.bin'];
    fid=fopen(ffnm);
    tmpdata=fread(fid,'ubit16');
    fclose(fid);
    fdr = reshape(tmpdata,c,r);
    slopepix = nan(c,r);
    if bid == 112% this basin has nan at the boundary, needs to be treated differently
        excep=4;
    else
        excep=2;
    end
    for i = 2:c-1
        for j = excep:r-1
%             if ismember(c*j+i,ind)
%                break 
%             end
            tmpz = fdr(i,j);
            tmp1 = i;tmp2 = j;
            d1 = dem(i,j);

            if tmpz == 1
                ik = tmp1+1;
                jk = tmp2-1;
                dis = 900*sqrt(2);
            elseif  tmpz == 2
                ik = tmp1+1;
                jk = tmp2;
                dis = 900;
            elseif  tmpz == 4
                ik = tmp1+1;
                jk = tmp2+1;
                dis = 900*sqrt(2);
            elseif  tmpz == 8
                ik = tmp1;
                jk = tmp2+1;
                dis = 900;
            elseif  tmpz == 16
                ik = tmp1-1;
                jk = tmp2+1;
                dis = 900*sqrt(2);
            elseif  tmpz == 32
                ik = tmp1-1;
                jk = tmp2;
                dis = 900;
            elseif  tmpz == 64
                ik = tmp1-1;
                jk = tmp2-1;
                dis = 900*sqrt(2);
            elseif  tmpz == 128
                ik = tmp1;
                jk = tmp2-1;
                dis = 900;
            end
            d2 = dem(ik,jk);
            slopepix(i,j) = (d1-d2)/dis;
            if slopepix(i,j)==0
                slopepix(i,j)=0.0001;
            end

        end
    end
    slopesum{bid,1} = slopepix;
    slopemeansum(bid,1) = mean(slopepix(intt));
end
%
dssum = {};
dosum = {};
pmeans = {};
pchanges = {};pospchanges = {};negpchanges = {};
monthnum = [];
yearnum = [];

for bid = [1:129 132:133 466:480 482:504 507:522]
    gid = WORLD_events(bid,1);
    tmp = load(['/shared/dondo/home/ml423/DA/Critical matrices/World_Basins/inttind/Basin' num2str(gid) '_ind.mat'],'ind');
    ind = tmp.ind;
    tmp = load(['/shared/dondo/home/ml423/DA/Critical matrices/World_Basins/inttind/Basin' num2str(gid) '_intt.mat'],'intt');
    intt = tmp.intt;       

    d_s = load(['/shared/dondo/home/ml423/DA/Critical matrices/World_Basins/distance/Basin' num2str(gid) 'dis_stream.mat']);
    d_s = d_s.dis_stream;
    d_o = load(['/shared/dondo/home/ml423/DA/Critical matrices/World_Basins/distance/Basin' num2str(gid) 'dis_outlet.mat']);
    d_o = d_o.dis_outlet;
    d_s(ind) = nan;d_o(ind) = nan;
    d_s = d_s+1;%dis_stream is set to be 1 for minimum
    d_s = d_s/max(d_s(:));
    d_o = d_o/max(d_o(:));
    
    for eid = 1:149
        if isnan(eqcontrol(bid,eid))||isempty(rain00_sum{bid,eid})||isempty(rain34_sum{bid,eid})
            continue     
        end
        tmp = rain00_sum{bid,eid};
        tmp34 = rain34_sum{bid,eid};
        tmpdsum = [];
        tmpdoum = [];
        pmeansum = [];
        pchangesum = [];
        pospchangesum = [];
        negpchangesum = [];
        for i = 1:size(tmp,3)
           tmps = tmp(:,:,i);
           tmps(ind) = nan;
           numeri = tmps.*d_s;
           tmpds = (sum(numeri(intt)))/sum(tmps(intt));
           tmpdsum = [tmpdsum,tmpds];
           
           numeri = tmps.*d_o;
           tmpdo = (sum(numeri(intt)))/sum(tmps(intt));
           tmpdoum = [tmpdoum,tmpdo];
           
           pmean = mean(tmps(intt));
           pmeansum = [pmeansum,pmean];
           
           tmps34 = tmp34(:,:,i);
           tmpdif = tmps34-tmps;
           
           pchange = mean(tmpdif(intt));
           pchangesum = [pchangesum,pchange];
           
           positivedif = find(tmpdif>0);
           negativedif = find(tmpdif<0);
           recon = zeros(size(tmps34,1),size(tmps34,2));
           recon(ind) = nan;
           recon(positivedif) = tmpdif(positivedif);
           numeri = recon.*d_s;
           tmpds = (sum(numeri(intt)))/sum(recon(intt));
           pospchangesum = [pospchangesum,tmpds];
           recon = zeros(size(tmps34,1),size(tmps34,2));
           recon(ind) = nan;
           recon(negativedif) = -1*tmpdif(negativedif);
           numeri = recon.*d_s;
           tmpds = (sum(numeri(intt)))/sum(recon(intt));
           negpchangesum = [negpchangesum,tmpds];
        end
        dssum{bid,eid} = tmpdsum;
        dosum{bid,eid} = tmpdoum;
        pmeans{bid,eid} = pmeansum;
        pchanges{bid,eid} = pchangesum;
        monthnum(bid,eid) = mod(floor(WORLD_events(bid,eid+1)/100),100);
        yearnum(bid,eid) = floor(WORLD_events(bid,eid+1)/10000);
        pospchanges{bid,eid} = pospchangesum;
        negpchanges{bid,eid} = negpchangesum;
    end
end

% d_area slopemeansum and the above three variables

% contruct test and all data for alps, tests ~ f(basin,event)


allalp = [1:129 132:133 466:480 482:504 507:522];
alleve = [];
% all alps events are below:

for bid = [1:129 132:133 466:480 482:504 507:522]
    if ~ismember(bid,problem_basins)
        for eid = 1:149
            if isnan(eqcontrol(bid,eid))||isempty(rain00_sum{bid,eid})||isempty(rain34_sum{bid,eid})
                continue
            end
            alleve = [alleve;[bid,eid]];
        end
    end
end
alleveori = alleve;

% the AI model correspond to kfold=1 is used for other mountains
for kfold = 1:5
    rng(42)
    if kfold<=4
        randloc = randperm(length(alleve),round(0.2*(length(alleveori))));
        basin_event_test = alleve(randloc,:);
        %save(['/shared/dondo/home/ml423/world_mts_runs/AI_test_basin_event.mat'],'basin_event_test')
        %save(['/shared/dondo/home/ml423/world_mts_runs/AI_test_basin_event_fold' num2str(kfold) '.mat'],'basin_event_test')
        alleve = setdiff(alleve,basin_event_test,'rows');
    else
        basin_event_test = alleve;
        %save(['/shared/dondo/home/ml423/world_mts_runs/AI_test_basin_event_fold' num2str(kfold) '.mat'],'basin_event_test')
    end


    dsall = [];dsall_test = [];
    doall = [];doall_test = [];
    pall = [];pall_test = [];
    demall = [];demall_test = [];
    sall = [];sall_test = [];mall = [];mall_test = [];yall = [];yall_test = [];
    aall = [];aall_test = [];
    cvfloodall = [];cvfloodall_test = [];
    floodmeanall = [];floodmeanall_test = [];floodmaxall = [];floodmaxall_test = [];
    floodall = [];floodall_test = [];
    targetall = [];targetall_test = [];
    pospdsall = [];pospds_test = [];
    negpdsall = [];negpds_test = [];

    for bid = [1:129 132:133 466:480 482:504 507:522]

        for eid = 1:149
            if ismember([bid,eid],basin_event_test,'rows')%&&~isempty(dssum{bid,eid})
                numhours = length(dssum{bid,eid});

                dsall_test = [dsall_test;dssum{bid,eid}'];
                doall_test = [doall_test;dosum{bid,eid}'];
                pall_test = [pall_test;pmeans{bid,eid}'];
                sall_test = [sall_test;slopemeansum(bid)*ones(numhours,1)];
                demall_test = [demall_test;demsum(bid)*ones(numhours,1)];
                mall_test = [mall_test;monthnum(bid,eid)*ones(numhours,1)];
                yall_test = [yall_test;yearnum(bid,eid)*ones(numhours,1)];
                aall_test = [aall_test;d_area(bid)*ones(numhours,1)];
                cvfloodall_test = [cvfloodall_test;cvflood(bid,eid)*ones(numhours,1)];
                floodmeanall_test = [floodmeanall_test;meanflood(bid,eid)*ones(numhours,1)];
                floodmaxall_test = [floodmaxall_test;maxflood(bid,eid)*ones(numhours,1)];
                floodall_test = [floodall_test;simu_flood{bid,eid}'];
                targetall_test = [targetall_test;pchanges{bid,eid}'];
                pospds_test = [pospds_test;pospchanges{bid,eid}'];
                negpds_test = [negpds_test;negpchanges{bid,eid}'];

            elseif ismember([bid,eid],alleveori,'rows')
                numhours = length(dssum{bid,eid});

                dsall = [dsall;dssum{bid,eid}'];
                doall = [doall;dosum{bid,eid}'];
                pall = [pall;pmeans{bid,eid}'];
                sall = [sall;slopemeansum(bid)*ones(numhours,1)];
                demall = [demall;demsum(bid)*ones(numhours,1)];
                mall = [mall;monthnum(bid,eid)*ones(numhours,1)];
                yall = [yall;yearnum(bid,eid)*ones(numhours,1)];
                aall = [aall;d_area(bid)*ones(numhours,1)];
                cvfloodall = [cvfloodall;cvflood(bid,eid)*ones(numhours,1)];
                floodmeanall = [floodmeanall;meanflood(bid,eid)*ones(numhours,1)];
                floodmaxall = [floodmaxall;maxflood(bid,eid)*ones(numhours,1)];
                floodall = [floodall;simu_flood{bid,eid}'];
                targetall = [targetall;pchanges{bid,eid}'];
                pospdsall = [pospdsall;pospchanges{bid,eid}'];
                negpdsall = [negpdsall;negpchanges{bid,eid}'];
            end

        end
    end





    figure
    scatter(pall,targetall,12,'filled')

    figure
    scatter(dsall,targetall,12,'filled')



    rng(38)
    % allin = [dsall';doall';pall';sall';aall';mall';yall';cvfloodall';floodmaxall';floodmeanall';floodall'];
    % targ = targetall';
    % net = feedforwardnet([12]);
    % net = train(net,allin,targ);
    % allin_test = [dsall_test';doall_test';pall_test';sall_test';aall_test';mall_test';yall_test';cvfloodall_test';floodmaxall_test';floodmeanall_test';floodall_test'];
    % predic1 = net(allin_test);% example case
    %

    % short version target delta P
    allin = [dsall';doall';pall';sall';aall';demall';cvfloodall';floodmaxall';floodmeanall';floodall'];
    targ = targetall';


    rng(39)
    net = feedforwardnet([12 12]);
    net = train(net,allin,targ);
    allin_test = [dsall_test';doall_test';pall_test';sall_test';aall_test';demall_test';cvfloodall_test';floodmaxall_test';floodmeanall_test';floodall_test'];
    predicnn = net(allin_test);% example case


    figure
    scatter(targetall_test,predicnn,14,'filled')

    %{
    % 20 ensemble AI models of Neural networks follow above, use basin P as a mean
    nn = [16];
    nn = [12 12];
    rng(37)
    net = feedforwardnet(nn);
    net1 = train(net,allin,targ);
    %predic1 = net1(allin_test);% example case
    rng(36)
    net = feedforwardnet(nn);
    net2 = train(net,allin,targ);
    %predic2 = net2(allin_test);% example case
    rng(35)
    net = feedforwardnet(nn);
    net3 = train(net,allin,targ);
    %predic3 = net3(allin_test);% example case
    rng(34)
    net = feedforwardnet(nn);
    net4 = train(net,allin,targ);
    %predic4 = net4(allin_test);% example case
    rng(33)
    net = feedforwardnet(nn);
    net5 = train(net,allin,targ);
    %predic5 = net5(allin_test);% example case
    rng(32)
    net = feedforwardnet(nn);
    net6 = train(net,allin,targ);
    %predic6 = net6(allin_test);% example case
    rng(31)
    net = feedforwardnet(nn);
    net7 = train(net,allin,targ);
    %predic7 = net7(allin_test);% example case
    rng(30)
    net = feedforwardnet(nn);
    net8 = train(net,allin,targ);
    %predic8 = net8(allin_test);% example case
    rng(29)
    net = feedforwardnet(nn);
    net9 = train(net,allin,targ);
    %predic9 = net9(allin_test);% example case
    rng(28)
    net = feedforwardnet(nn);
    net10 = train(net,allin,targ);
    %predic10 = net10(allin_test);% example case
    rng(27)
    net = feedforwardnet(nn);
    net11 = train(net,allin,targ);
    %predic11 = net11(allin_test);% example case
    rng(26)
    net = feedforwardnet(nn);
    net12 = train(net,allin,targ);
    %predic12 = net12(allin_test);% example case
    rng(25)
    net = feedforwardnet(nn);
    net13 = train(net,allin,targ);
    %predic13 = net13(allin_test);% example case
    rng(24)
    net = feedforwardnet(nn);
    net14 = train(net,allin,targ);
    %predic14 = net14(allin_test);% example case
    rng(23)
    net = feedforwardnet(nn);
    net15 = train(net,allin,targ);
    %predic15 = net15(allin_test);% example case
    rng(22)
    net = feedforwardnet(nn);
    net16 = train(net,allin,targ);
    %predic16 = net16(allin_test);% example case
    rng(21)
    net = feedforwardnet(nn);
    net17 = train(net,allin,targ);
    %predic17 = net17(allin_test);% example case
    rng(20)
    net = feedforwardnet(nn);
    net18 = train(net,allin,targ);
    %predic18 = net18(allin_test);% example case
    rng(19)
    net = feedforwardnet(nn);
    net19 = train(net,allin,targ);
    %predic19 = net19(allin_test);% example case
    rng(18)
    net = feedforwardnet(nn);
    net20 = train(net,allin,targ);
    %predic20 = net20(allin_test);% example case

    save(['/shared/dondo/home/ml423/world_mts_runs/AI_test_basin_event_fold_' num2str(kfold) '_net' num2str(1,'%2.2d') '.mat'],'net1')
    save(['/shared/dondo/home/ml423/world_mts_runs/AI_test_basin_event_fold_' num2str(kfold) '_net' num2str(2,'%2.2d') '.mat'],'net2')
    save(['/shared/dondo/home/ml423/world_mts_runs/AI_test_basin_event_fold_' num2str(kfold) '_net' num2str(3,'%2.2d') '.mat'],'net3')
    save(['/shared/dondo/home/ml423/world_mts_runs/AI_test_basin_event_fold_' num2str(kfold) '_net' num2str(4,'%2.2d') '.mat'],'net4')
    save(['/shared/dondo/home/ml423/world_mts_runs/AI_test_basin_event_fold_' num2str(kfold) '_net' num2str(5,'%2.2d') '.mat'],'net5')
    save(['/shared/dondo/home/ml423/world_mts_runs/AI_test_basin_event_fold_' num2str(kfold) '_net' num2str(6,'%2.2d') '.mat'],'net6')
    save(['/shared/dondo/home/ml423/world_mts_runs/AI_test_basin_event_fold_' num2str(kfold) '_net' num2str(7,'%2.2d') '.mat'],'net7')
    save(['/shared/dondo/home/ml423/world_mts_runs/AI_test_basin_event_fold_' num2str(kfold) '_net' num2str(8,'%2.2d') '.mat'],'net8')
    save(['/shared/dondo/home/ml423/world_mts_runs/AI_test_basin_event_fold_' num2str(kfold) '_net' num2str(9,'%2.2d') '.mat'],'net9')
    save(['/shared/dondo/home/ml423/world_mts_runs/AI_test_basin_event_fold_' num2str(kfold) '_net' num2str(10,'%2.2d') '.mat'],'net10')
    save(['/shared/dondo/home/ml423/world_mts_runs/AI_test_basin_event_fold_' num2str(kfold) '_net' num2str(11,'%2.2d') '.mat'],'net11')
    save(['/shared/dondo/home/ml423/world_mts_runs/AI_test_basin_event_fold_' num2str(kfold) '_net' num2str(12,'%2.2d') '.mat'],'net12')
    save(['/shared/dondo/home/ml423/world_mts_runs/AI_test_basin_event_fold_' num2str(kfold) '_net' num2str(13,'%2.2d') '.mat'],'net13')
    save(['/shared/dondo/home/ml423/world_mts_runs/AI_test_basin_event_fold_' num2str(kfold) '_net' num2str(14,'%2.2d') '.mat'],'net14')
    save(['/shared/dondo/home/ml423/world_mts_runs/AI_test_basin_event_fold_' num2str(kfold) '_net' num2str(15,'%2.2d') '.mat'],'net15')
    save(['/shared/dondo/home/ml423/world_mts_runs/AI_test_basin_event_fold_' num2str(kfold) '_net' num2str(16,'%2.2d') '.mat'],'net16')
    save(['/shared/dondo/home/ml423/world_mts_runs/AI_test_basin_event_fold_' num2str(kfold) '_net' num2str(17,'%2.2d') '.mat'],'net17')
    save(['/shared/dondo/home/ml423/world_mts_runs/AI_test_basin_event_fold_' num2str(kfold) '_net' num2str(18,'%2.2d') '.mat'],'net18')
    save(['/shared/dondo/home/ml423/world_mts_runs/AI_test_basin_event_fold_' num2str(kfold) '_net' num2str(19,'%2.2d') '.mat'],'net19')
    save(['/shared/dondo/home/ml423/world_mts_runs/AI_test_basin_event_fold_' num2str(kfold) '_net' num2str(20,'%2.2d') '.mat'],'net20')

    %}
    %{
    % apply NN to cross-validation alps
    predic20 = net20(allin_test);% example case
    predic19 = net19(allin_test);% example case
    predic18 = net18(allin_test);% example case
    predic17 = net17(allin_test);% example case
    predic16 = net16(allin_test);% example case
    predic15 = net15(allin_test);% example case
    predic14 = net14(allin_test);% example case
    predic13 = net13(allin_test);% example case
    predic12 = net12(allin_test);% example case
    predic11 = net11(allin_test);% example case
    predic10 = net10(allin_test);% example case
    predic9 = net9(allin_test);% example case
    predic8 = net8(allin_test);% example case
    predic7 = net7(allin_test);% example case
    predic6 = net6(allin_test);% example case
    predic5 = net5(allin_test);% example case
    predic4 = net4(allin_test);% example case
    predic3 = net3(allin_test);% example case
    predic2 = net2(allin_test);% example case
    predic1 = net1(allin_test);% example case
    
    
    % apply the AI model to alps f(basin,event)
    maxrep = 20;

    accuhour = 0;
    for bid = [1:129 132:133 466:480 482:504 507:522]
        gid = WORLD_events(bid,1);
        gloc = find(Basins_sum(:,1)==gid);

        r = Basins_sum(gloc,6);
        c = Basins_sum(gloc,7);

        for eid = 1:149
            if ismember([bid,eid],basin_event_test,'rows')
                event = WORLD_events(bid,eid+1);

                ffnm=['/shared/dondo/home/ml423/world_mts/Precip/Event/Basin' num2str(gid) '_input_ERA5_900m_FIX/' num2str(event) '/Precipitation_ERA5_land_oriunit'];
                if exist(ffnm)
                    orir = ['/shared/dondo/home/ml423/world_mts/Precip/Event/Basin' num2str(gid) '_input_ERA5_900m_FIX/' num2str(event) '/'];
                    out_dir2 = ['/shared/dondo/home/ml423/world_mts/Basin' num2str(gid) '_input_ERA5_900m_tmp/' num2str(event) '/'];
                    copyfile([orir '*'],out_dir2)
                    disp('transfer done')
                else
                    disp('error')
                end

                ffnm=['/shared/dondo/home/ml423/world_mts/Basin' num2str(gid) '_input_ERA5_900m_tmp/' num2str(event) '/Precipitation_ERA5_land_oriunit'];
                fid=fopen(ffnm,'rb','ieee-le');
                tmpdata=fread(fid,inf,'single');
                fclose(fid);
                NODAc2o5=reshape(tmpdata,c,r,720);

                kind = eval_ind{bid,eid};

                for replic = 1:maxrep
                    predics = eval(['predic' num2str(replic)]);
                    %predics = targetall_test;% actually try this
                    reviserain = [];


                    for k = 1:720
                        tmpori = NODAc2o5(:,:,k);
                        if ismember(k,kind)
                            accuhour = accuhour+1;
                            if replic>1&&k==kind(1)
                                accuhour = accuhour-length(kind);
                            end
                            tmpupdate = tmpori+predics(accuhour);
                            reviserain(:,:,k) = tmpupdate;
                        else
                            reviserain(:,:,k) = tmpori;
                        end
                        %[k,replic,accuhour]
                    end
                    %bid
                    %eid
                    accuhour
                    out_dirp = ['/shared/dondo/home/ml423/world_mts/Basin' num2str(gid) '_input_ERA5_900m_tmp/' num2str(event) '/'];
                    mkdir(out_dirp)

                    reviserain(reviserain<0) = 0;
                    fnm = [out_dirp 'Precipitation_ERA5_land_ori_' num2str(replic,'%3.3d')];
                    raindata = reshape(reviserain,[c,r,720]);
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
            end
        end
    end
    %}
    
    % build Random forest regressor below
    X_reg = allin';X_reg = X_reg(:,[1:5 7 10]);
    Y_reg = targ';
    X_test = allin_test';X_test = X_test(:,[1:5 7 10]);

    disp('--- Regression Example ---');
    rng('default'); % For reproducibility
    numTrees_reg = 30;
    Mdl_reg1 = TreeBagger(numTrees_reg, X_reg, Y_reg, ...
        'Method', 'regression', ...
        'OOBPrediction', 'on', ...
        'OOBPredictorImportance', 'on', ...
        'MinLeafSize', 10); 
    rng(45); % For reproducibility, 0 42 43 44 45
    % Train a TreeBagger model for regression
    numTrees_reg = 30;
    Mdl_reg5 = TreeBagger(numTrees_reg, X_reg, Y_reg, ...
        'Method', 'regression', ...
        'OOBPrediction', 'on', ...
        'OOBPredictorImportance', 'on', ...
        'MinLeafSize', 10); % Enable out-of-bag error estimation
    rng(44); % For reproducibility, 0 42 43 44 45
    numTrees_reg = 30;
    Mdl_reg4 = TreeBagger(numTrees_reg, X_reg, Y_reg, ...
        'Method', 'regression', ...
        'OOBPrediction', 'on', ...
        'OOBPredictorImportance', 'on', ...
        'MinLeafSize', 10); % Enable out-of-bag error estimation
    rng(43); % For reproducibility, 0 42 43 44 45
    numTrees_reg = 30;
    Mdl_reg3 = TreeBagger(numTrees_reg, X_reg, Y_reg, ...
        'Method', 'regression', ...
        'OOBPrediction', 'on', ...
        'OOBPredictorImportance', 'on', ...
        'MinLeafSize', 10); % Enable out-of-bag error estimation
    rng(42); % For reproducibility, 0 42 43 44 45
    numTrees_reg = 30;
    Mdl_reg2 = TreeBagger(numTrees_reg, X_reg, Y_reg, ...
        'Method', 'regression', ...
        'OOBPrediction', 'on', ...
        'OOBPredictorImportance', 'on', ...
        'MinLeafSize', 10); % Enable out-of-bag error estimation


    predic1alps = predict(Mdl_reg1, X_test);
    predic2alps = predict(Mdl_reg2, X_test);
    predic3alps = predict(Mdl_reg3, X_test);
    predic4alps = predict(Mdl_reg4, X_test);
    predic5alps = predict(Mdl_reg5, X_test);
    
    % apply the AI model to alps f(basin,event) for cross validation
    maxrep = 5;

    accuhour = 0;
    for bid = [1:129 132:133 466:480 482:504 507:522]
        gid = WORLD_events(bid,1);
        gloc = find(Basins_sum(:,1)==gid);

        r = Basins_sum(gloc,6);
        c = Basins_sum(gloc,7);

        for eid = 1:149
            if ismember([bid,eid],basin_event_test,'rows')
                event = WORLD_events(bid,eid+1);

                ffnm=['/shared/dondo/home/ml423/world_mts/Precip/Event/Basin' num2str(gid) '_input_ERA5_900m_FIX/' num2str(event) '/Precipitation_ERA5_land_oriunit'];
                if exist(ffnm)
                    orir = ['/shared/dondo/home/ml423/world_mts/Precip/Event/Basin' num2str(gid) '_input_ERA5_900m_FIX/' num2str(event) '/'];
                    out_dir2 = ['/shared/dondo/home/ml423/world_mts/Basin' num2str(gid) '_input_ERA5_900m_tmp/' num2str(event) '/'];
                    copyfile([orir '*'],out_dir2)
                    disp('transfer done')
                else
                    disp('error')
                end

                ffnm=['/shared/dondo/home/ml423/world_mts/Basin' num2str(gid) '_input_ERA5_900m_tmp/' num2str(event) '/Precipitation_ERA5_land_oriunit'];
                fid=fopen(ffnm,'rb','ieee-le');
                tmpdata=fread(fid,inf,'single');
                fclose(fid);
                NODAc2o5=reshape(tmpdata,c,r,720);

                kind = eval_ind{bid,eid};

                for replic = 1:maxrep
                    predics = eval(['predic' num2str(replic) 'alps']);
                    %predics = targetall_test;% actually try this
                    reviserain = [];


                    for k = 1:720
                        tmpori = NODAc2o5(:,:,k);
                        if ismember(k,kind)
                            accuhour = accuhour+1;
                            if replic>1&&k==kind(1)
                                accuhour = accuhour-length(kind);
                            end
                            tmpupdate = tmpori+predics(accuhour);
                            reviserain(:,:,k) = tmpupdate;
                        else
                            reviserain(:,:,k) = tmpori;
                        end
                        %[k,replic,accuhour]
                    end
                    %bid
                    %eid
                    accuhour
                    out_dirp = ['/shared/dondo/home/ml423/world_mts/Basin' num2str(gid) '_input_ERA5_900m_tmp/' num2str(event) '/'];
                    %mkdir(out_dirp)

                    reviserain(reviserain<0) = 0;
                    fnm = [out_dirp 'Precipitation_ERA5_lRFR_ori_' num2str(replic,'%3.3d')];
                    raindata = reshape(reviserain,[c,r,720]);
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
            end
        end
    end

end

%% Apply Alps model to basins in the Andes (cross-regional validation)
% calculate slopes
andessmall = find(gauge_events_first_last(:,5)<1200);
adsall = [];
for ki = 1:length(andessmall)
    bid = find(WORLD_events(:,1)==gauge_events_first_last(andessmall(ki),1));
    if bid>=134&bid<=465
        if gauge_events_first_last(andessmall(ki),3)>19900000
            adsall = [adsall;bid];
        end
    end
end
length(adsall)


% slopemeansumandes = [];
% demsumandes = [];
for andis = 1:length(adsall)
    bid = adsall(andis);
    gid = WORLD_events(bid,1);
    gloc = find(Basins_sum(:,1)==gid);
    
    r = Basins_sum(gloc,6);
    c = Basins_sum(gloc,7);
    tmp = load(['/shared/dondo/home/ml423/DA/Critical matrices/World_Basins/inttind/Basin' num2str(gid) '_ind.mat'],'ind');
    ind = tmp.ind;
    tmp = load(['/shared/dondo/home/ml423/DA/Critical matrices/World_Basins/inttind/Basin' num2str(gid) '_intt.mat'],'intt');
    intt = tmp.intt; 
    ffnm=['/shared/dondo/home/ml423/world_mts/Basin' num2str(gid) '_fixed_900m/dem.bin'];
    fid=fopen(ffnm,'rb','ieee-le');
    tmpdata=fread(fid,inf,'single');
    fclose(fid);
    dem = reshape(tmpdata,c,r);
    %dem(ind) = nan;
    demsumandes(bid) = mean(dem(intt));

    ffnm=['/shared/dondo/home/ml423/world_mts/Basin' num2str(gid) '_fixed_900m/fdr.bin'];
    fid=fopen(ffnm);
    tmpdata=fread(fid,'ubit16');
    fclose(fid);
    fdr = reshape(tmpdata,c,r);
    slopepix = nan(c,r);
    if bid == 341
        excep=3;
    else
        excep=2;
    end
    for i = excep:c-1
        for j = 2:r-1
%             if ismember(c*j+i,ind)
%                break 
%             end
            tmpz = fdr(i,j);
            tmp1 = i;tmp2 = j;
            d1 = dem(i,j);

            if tmpz == 1
                ik = tmp1+1;
                jk = tmp2-1;
                dis = 900*sqrt(2);
            elseif  tmpz == 2
                ik = tmp1+1;
                jk = tmp2;
                dis = 900;
            elseif  tmpz == 4
                ik = tmp1+1;
                jk = tmp2+1;
                dis = 900*sqrt(2);
            elseif  tmpz == 8
                ik = tmp1;
                jk = tmp2+1;
                dis = 900;
            elseif  tmpz == 16
                ik = tmp1-1;
                jk = tmp2+1;
                dis = 900*sqrt(2);
            elseif  tmpz == 32
                ik = tmp1-1;
                jk = tmp2;
                dis = 900;
            elseif  tmpz == 64
                ik = tmp1-1;
                jk = tmp2-1;
                dis = 900*sqrt(2);
            elseif  tmpz == 128
                ik = tmp1;
                jk = tmp2-1;
                dis = 900;
            end
            d2 = dem(ik,jk);
            slopepix(i,j) = (d1-d2)/dis;
            if slopepix(i,j)==0
                slopepix(i,j)=0.0001;
            end

        end
    end
    slopemeansumandes(bid,1) = mean(slopepix(intt));
end
%


flowtype = {'streamflowarea','interflowarea','overlandflowarea','baseflowarea'};
ftp = 1;

% kge_opt = [];kge_ori = [];
% rain00_sum = {};
% rain34_sum = {};
% eval_ind = {};
%simu_floodandes = {};
%cvflood = nan(133,149);maxflood = nan(133,149);meanflood = nan(133,149);
for andis = 1:length(adsall)
    bid = adsall(andis);
    gid = WORLD_events(bid,1);
    
    gloc = find(Basins_sum(:,1)==gid);
    
    r = Basins_sum(gloc,6);
    c = Basins_sum(gloc,7);
    rg = Basins_sum(gloc,8);
    cg = Basins_sum(gloc,9);
    d_areaandes(bid) = Basins_sum(gloc,5);
    
    tmp = load(['/shared/dondo/home/ml423/DA/Critical matrices/World_Basins/inttind/Basin' num2str(gid) '_ind.mat'],'ind');
    ind = tmp.ind;
    tmp = load(['/shared/dondo/home/ml423/DA/Critical matrices/World_Basins/inttind/Basin' num2str(gid) '_intt.mat'],'intt');
    intt = tmp.intt;
    ffnm=['/shared/dondo/home/ml423/world_mts/Basin' num2str(gid) '_fixed_900m/dem.bin'];
    fid=fopen(ffnm,'rb','ieee-le');
    tmpdata=fread(fid,inf,'single');
    fclose(fid);
    dem = reshape(tmpdata,c,r);
    dem(ind) = nan;
    

    ICn = 4140;
    ICC = 0;% IF ICC is part of the algo
    
    fileind = -100000*ICn;

    
    ixx = find(~isnan(WORLD_events(bid,:)));
    etot = length(ixx)-1;
    if etot == 0
        continue
    else
        
        for eid = 1:etot
            close all
            if isnan(event_quality(bid,eid+1))
                disp('Event quality is bad, peak too early or too late')
                continue
            else
                rainevent = num2str(WORLD_events(bid,eid+1));
                yth = floor(WORLD_events(bid,eid+1)/10000);
                dth = mod(WORLD_events(bid,eid+1),100);
                mth = floor((WORLD_events(bid,eid+1)-yth*10000)/100);
                
                thresdate = datetime(yth,mth,dth);
                
                durati = (thresdate-days(15)):days(1):(thresdate+days(15));
                [dury,durm,durd] = ymd(durati);
                tst1 = num2str(durm(1)*100+durd(1),'%4.4d');
                tst2 = num2str(durm(30)*100+durd(30),'%4.4d');
                
                % the condition below has to change, most of events in
                % Andes is not used IRCICC because no computational resouce
                ffnm=['/shared/dondo/home/ml423/world_mts_runs/IRC1/Basin' num2str(gid) 'outputs/refxtmp' num2str(str2num(rainevent)*100+fileind+0) 'fullresults/' rainevent '_SW/' rainevent '/streamflowarea'];
                fid=fopen(ffnm,'rb','ieee-le');
                if fid == -1
                    continue
                else
                    fclose(fid);
                    
                    ffnm = ['/shared/dondo/home/ml423/world_mts_runs/Basins_events_chars/'];
                    fid = fopen([ffnm 'Basin' num2str(gid) '_event' rainevent '.txt'],'r');
                    if fid == -1
                        fclose(fid);
                        continue
                    end
                    windowtimes = fscanf(fid,'%d\n',[4,1]);
                    fclose(fid);
                    if length(windowtimes)<4
                        continue 
                    end
                    %close all
       
                    
                    if 1==1
                        
                        %fh = APL_frontcut(bid,eid);
                        
                        tstr = [num2str(dury(1)) '-' num2str(durm(1),'%2.2d') '/' num2str(durd(1),'%2.2d') '-' num2str(durm(30),'%2.2d') '/' num2str(durd(30),'%2.2d')];
                        
                        
                        obstep = 90;
                        
                        tmp = load(['/shared/dondo/home/ml423/world_mts/obs_events/events/obs_' num2str(gid) '_event' rainevent '.mat']);
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
                        
                        
                        
                        tmp = load(['/shared/dondo/home/ml423/world_mts/Precip/Event/ERA5_Basin' num2str(gid) '_' rainevent '.mat']);
                        s4dbkc = tmp.rain;
                        
                        MCp = zeros(c,r);
                        pts = [];
                        dbkcr = 0;rmax = [];
                        for i = 1:720
                            MCp = MCp + s4dbkc{i};
                            pts(i) = mean(s4dbkc{i}(intt));
                            dbkcr = dbkcr+pts(i);
                            rmax(i) = max(s4dbkc{i}(intt));
                        end
                        MCp(ind) = nan;
                        ori_mean = nanmean(MCp(:));
                        ama = jet(40);ama(1,:) = [1 1 1];
               
                        
                        % for Nature Science paper, consider the major events only
                        
                        launchp = find(rmax>0.5);% when rainfall >1mm/hr launch particles
                        ptmaxthre = 0.5*max(rmax);% used as a threshold to separate events
                        la2 = find(launchp<(windowtimes(4)*8));
                        la1 = find(launchp<(windowtimes(3)*8));
                        accu0 = 0;
                        if isempty(la1)
                            continue
                        end
                        for tmpla = launchp(la1(end)):-1:1
                            if ismember(tmpla,launchp)
                                accu0 = 0;
                                rainbegintime = tmpla;
                            else
                                accu0 = accu0+1;
                            end
                            
                            if area>1000
                                noraininterval = 12;% for big basins, it is likely there is always rain in the basin, therefore using a small norain limit to distinguish different rain events.
                            else
                                noraininterval = 24;
                            end
                            if accu0>noraininterval
                                % consecutive x hours of basin max rain less than 0.5mm/hr,
                                % this should be counted as a different rain cell, should not
                                % be tracked!
                                % one exception is if another major rain is 48-24 ahead, should be
                                % considered as well
                                rsmallest = max([(tmpla-noraininterval),1]);

                                reconsider_period = rsmallest:tmpla;
                                %reconsider_period = (tmpla-noraininterval):tmpla;
                                secondary_rain = find(rmax(reconsider_period)>0.5);
                                
                                if max(rmax(reconsider_period))>ptmaxthre
                                    % adjust the rainbegintime here
                                    rainbegintime = reconsider_period(secondary_rain(1));
                                end
                                smalltmpla = max([floor(tmpla/8),2]);
                                if obsnow(smalltmpla)>obsnow(smalltmpla-1)
                                    % meaning obsflow is still rising.
                                    if ~isempty(secondary_rain)
                                        rainbegintime = reconsider_period(secondary_rain(1));
                                    end
                                end
                                
                                break
                            end
                        end
                        
                        la0 = find(launchp == rainbegintime);
                        
                        refinedlaunchp = launchp(la0:la2(end));
                        stat_left_bd_hydro = windowtimes(1);% systematic defined
                        stat_left_bd_rain = floor(rainbegintime/8);%target rainevent adjusted
                        stat_left_bd = max([stat_left_bd_hydro,stat_left_bd_rain]);

                        upcap = min([90 (windowtimes(4)+6)]);

                        if upcap<=stat_left_bd
                            continue
                        end

                        shj = [1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,6,6,6,7,7,7];
                        tww = 3;
                        KGEwxy = [];
                        no3sum = [];
                        ftot = [];
                        for p = 0
                            
                            %itranum = 0;
                            
                            if p==0
                                filenam = str2num(rainevent)*100+fileind;
                                wnum = mod(filenam,10);
                                itranum = 0;
                            else
                                wnum = mod(filenam,10);
                                itranum = shj(p);
                                
                                if mod(p,tww)==1%window 2
                                    zpo = windowtimes(1);
                                elseif mod(p,tww)==2%window 3
                                    zpo = windowtimes(2);
                                elseif mod(p,tww)==0%window 4
                                    zpo = windowtimes(3);
                                end
                                
                            end
 
                            
                            ffnmx=['/shared/dondo/home/ml423/world_mts_runs/IRC1/Basin' num2str(gid) 'outputs/refxtmp' num2str(filenam) 'fullresults/' rainevent '_SW/' rainevent '/' flowtype{ftp}];
                            fid=fopen(ffnmx,'rb','ieee-le');
                            
                            if fid==-1
                                %    realIRCstat = [realIRCstat;[raininfob1(rit,1),nan,nan,nan,nan,nan]];
                                %disp('error')
                                %ffnm
                                
                                % meaning the files are not there, so move
                                % on to next iteration. THis can be because
                                % the user delete intermediate files to
                                % save space.
                                if p==0
                                    filenam = filenam+2;
                                else
                                    if mod(p,tww)==1
                                        filenam = filenam+1;
                                    elseif mod(p,tww)==2
                                        filenam = filenam+1;
                                    elseif mod(p,tww)==0
                                        filenam = filenam+8;
                                    end
                                end

                                if p<9
                                    continue
                                else
                                    break
                                end
                            end
                            tmpdata=fread(fid,inf,'single');
                            fclose(fid);
                            NODAc2o5=reshape(tmpdata,c,r,720);
                            flownoda5 = [];
                            for k = 1:720
                                flownoda5(k) = abs(NODAc2o5(cg,rg,k));
                            end
                            
                            no3 = flownoda5;
                            no3sum = [no3sum;no3];
                            
                            if p==0
                                nostivdbkc = no3;
                            end
                            
                            
                            ffnm=['/shared/dondo/home/ml423/world_mts_runs/IRC1/Basin' num2str(gid) 'outputs/refxtmp' num2str(filenam) 'fullresults/' rainevent '_SW/' rainevent '/precipitation'];
                            fid=fopen(ffnm,'rb','ieee-le');
                            tmpdata=fread(fid,inf,'single');
                            fclose(fid);
                            NODAc2o5=reshape(tmpdata,c,r,720);
                            flownoda5 = [];
                            meanr250 = 0;
                            tmpsum = zeros(c,r);
                            for k = 1:720
                                tmp = NODAc2o5(:,:,k);
                                tmpsum = tmpsum+tmp;
                                flownoda5(k) = mean(tmp(intt));
                                meanr250 = meanr250+flownoda5(k);
                            end
                            no2 = flownoda5;%GREEN
                            rainsd = std(tmpsum(intt));
                            
                            cal_window = (stat_left_bd-1)*8+1:(upcap-1)*8+1;
                            
                            if p==0
                                tstmp = find(flownoda5>0.1);% only consider rain>0.1mm/h in the original field
                                % ircicc rain>0.1 hours are a little bit
                                % more than original rain>0.1 hours have to
                                % be consistent for the number of hours for
                                % error modeling
                                rts1 = intersect(tstmp,cal_window);
                                rain00_sumandes{bid,eid} = NODAc2o5(:,:,rts1);
                                eval_indandes{bid,eid} = rts1;
                                simu_floodandes{bid,eid} = no3(rts1);
                                if isempty(rts1)
                                    cvfloodandes(bid,eid) = nan;
                                    maxfloodandes(bid,eid) = nan;
                                    meanfloodandes(bid,eid) = nan;
                                else
                                    cvfloodandes(bid,eid) = mean(no3(rts1))/std(no3(rts1));
                                    maxfloodandes(bid,eid) = max(no3(rts1));
                                    meanfloodandes(bid,eid) = mean(no3(rts1));
                                end
                            end
                            %
       
                            ffnm=['/shared/dondo/home/ml423/world_mts_runs/IRC1/Basin' num2str(gid) 'outputs/refxtmp' num2str(filenam+24) 'fullresults/' rainevent '_SW/' rainevent '/precipitation'];
                            fid=fopen(ffnm,'rb','ieee-le');
                            if fid == -1
                                rain34_sumandes{bid,eid} = nan(c,r,length(rts1));
                            else
                                
                            tmpdata=fread(fid,inf,'single');
                            fclose(fid);
                            NODAc2o5=reshape(tmpdata,c,r,720);
                            rain34_sumandes{bid,eid} = NODAc2o5(:,:,rts1);
                            end

                        end
                  
                    end
                end
            end
            [bid,eid]
        end
        
    end
    bid
end

%% Apply Alps models to Andes events that are created in another data server /taiga/ML/AI_20/Basinxoutputs/ due to data storage limit
% SQ: using ERA5L original as inputs, the outputs are in /taiga.

andessmall = find(gauge_events_first_last(:,5)<1200);
adsall = [];
for ki = 1:length(andessmall)
    bid = find(WORLD_events(:,1)==gauge_events_first_last(andessmall(ki),1));
    if bid>=134&bid<=465
        if gauge_events_first_last(andessmall(ki),3)>19900000
            adsall = [adsall;bid];
        end
    end
end
length(adsall)
flowtype = {'streamflowarea','interflowarea','overlandflowarea','baseflowarea'};
ftp = 1;

% kge_opt = [];kge_ori = [];
% rain00_sum = {};
% rain34_sum = {};
% eval_ind = {};
% simu_floodandes = {};
% cvflood = nan(133,149);maxflood = nan(133,149);meanflood = nan(133,149);
for andis = 1:length(adsall)
    bid = adsall(andis);
    gid = WORLD_events(bid,1);
    
    gloc = find(Basins_sum(:,1)==gid);
    
    r = Basins_sum(gloc,6);
    c = Basins_sum(gloc,7);
    rg = Basins_sum(gloc,8);
    cg = Basins_sum(gloc,9);
    d_areaandes(bid) = Basins_sum(gloc,5);
    
    tmp = load(['/shared/dondo/home/ml423/DA/Critical matrices/World_Basins/inttind/Basin' num2str(gid) '_ind.mat'],'ind');
    ind = tmp.ind;
    tmp = load(['/shared/dondo/home/ml423/DA/Critical matrices/World_Basins/inttind/Basin' num2str(gid) '_intt.mat'],'intt');
    intt = tmp.intt;
    ffnm=['/shared/dondo/home/ml423/world_mts/Basin' num2str(gid) '_fixed_900m/dem.bin'];
    fid=fopen(ffnm,'rb','ieee-le');
    tmpdata=fread(fid,inf,'single');
    fclose(fid);
    dem = reshape(tmpdata,c,r);
    dem(ind) = nan;
    
    ixx = find(~isnan(WORLD_events(bid,:)));
    etot = length(ixx)-1;
    if etot == 0
        continue
    else
        
        for eid = 1:etot
            close all
            if isnan(event_quality(bid,eid+1))
                disp('Event quality is bad, peak too early or too late')
                continue
            else
                rainevent = num2str(WORLD_events(bid,eid+1));
                yth = floor(WORLD_events(bid,eid+1)/10000);
                dth = mod(WORLD_events(bid,eid+1),100);
                mth = floor((WORLD_events(bid,eid+1)-yth*10000)/100);
                
                thresdate = datetime(yth,mth,dth);
                
                durati = (thresdate-days(15)):days(1):(thresdate+days(15));
                [dury,durm,durd] = ymd(durati);
                tst1 = num2str(durm(1)*100+durd(1),'%4.4d');
                tst2 = num2str(durm(30)*100+durd(30),'%4.4d');
                
                % the condition below has to change, most of events in
                % Andes is not used IRCICC because no computational resouce
                ffnm=['/taiga/ML/AI_20/Basin' num2str(gid) 'outputs/' rainevent '_SQ/' rainevent '/streamflowarea'];
                fid=fopen(ffnm,'rb','ieee-le');
                if fid == -1
                    continue
                else
                    fclose(fid);
                    
                    ffnm = ['/shared/dondo/home/ml423/world_mts_runs/Basins_events_chars/'];
                    fid = fopen([ffnm 'Basin' num2str(gid) '_event' rainevent '.txt'],'r');
                    if fid == -1
                        fclose(fid);
                        continue
                    end
                    windowtimes = fscanf(fid,'%d\n',[4,1]);
                    fclose(fid);
                    if length(windowtimes)<4
                        continue
                    end
                    %close all
                    
                    
                    if 1==1
                        
                        %fh = APL_frontcut(bid,eid);
                        
                        tstr = [num2str(dury(1)) '-' num2str(durm(1),'%2.2d') '/' num2str(durd(1),'%2.2d') '-' num2str(durm(30),'%2.2d') '/' num2str(durd(30),'%2.2d')];
                        
                        
                        obstep = 90;
                        
                        tmp = load(['/shared/dondo/home/ml423/world_mts/obs_events/events/obs_' num2str(gid) '_event' rainevent '.mat']);
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
                        
                        
                        
                        tmp = load(['/shared/dondo/home/ml423/world_mts/Precip/Event/ERA5_Basin' num2str(gid) '_' rainevent '.mat']);
                        s4dbkc = tmp.rain;
                        
                        MCp = zeros(c,r);
                        pts = [];
                        dbkcr = 0;rmax = [];
                        for i = 1:720
                            MCp = MCp + s4dbkc{i};
                            pts(i) = mean(s4dbkc{i}(intt));
                            dbkcr = dbkcr+pts(i);
                            rmax(i) = max(s4dbkc{i}(intt));
                        end
                        MCp(ind) = nan;
                        ori_mean = nanmean(MCp(:));
                        ama = jet(40);ama(1,:) = [1 1 1];
                        
                        
                        % for Nature Science paper, consider the major events only
                        
                        launchp = find(rmax>0.5);% when rainfall >1mm/hr launch particles
                        ptmaxthre = 0.5*max(rmax);% used as a threshold to separate events
                        la2 = find(launchp<(windowtimes(4)*8));
                        la1 = find(launchp<(windowtimes(3)*8));
                        accu0 = 0;
                        if isempty(la1)
                            continue
                        end
                        for tmpla = launchp(la1(end)):-1:1
                            if ismember(tmpla,launchp)
                                accu0 = 0;
                                rainbegintime = tmpla;
                            else
                                accu0 = accu0+1;
                            end
                            
                            if area>1000
                                noraininterval = 12;% for big basins, it is likely there is always rain in the basin, therefore using a small norain limit to distinguish different rain events.
                            else
                                noraininterval = 24;
                            end
                            if accu0>noraininterval
                                % consecutive x hours of basin max rain less than 0.5mm/hr,
                                % this should be counted as a different rain cell, should not
                                % be tracked!
                                % one exception is if another major rain is 48-24 ahead, should be
                                % considered as well
                                rsmallest = max([(tmpla-noraininterval),1]);
                                
                                reconsider_period = rsmallest:tmpla;
                                %reconsider_period = (tmpla-noraininterval):tmpla;
                                secondary_rain = find(rmax(reconsider_period)>0.5);
                                
                                if max(rmax(reconsider_period))>ptmaxthre
                                    % adjust the rainbegintime here
                                    rainbegintime = reconsider_period(secondary_rain(1));
                                end
                                smalltmpla = max([floor(tmpla/8),2]);
                                if obsnow(smalltmpla)>obsnow(smalltmpla-1)
                                    % meaning obsflow is still rising.
                                    if ~isempty(secondary_rain)
                                        rainbegintime = reconsider_period(secondary_rain(1));
                                    end
                                end
                                
                                break
                            end
                        end
                        
                        la0 = find(launchp == rainbegintime);
                        
                        refinedlaunchp = launchp(la0:la2(end));
                        stat_left_bd_hydro = windowtimes(1);% systematic defined
                        stat_left_bd_rain = floor(rainbegintime/8);%target rainevent adjusted
                        stat_left_bd = max([stat_left_bd_hydro,stat_left_bd_rain]);
                        
                        upcap = min([90 (windowtimes(4)+6)]);
                        
                        if upcap<=stat_left_bd
                            continue
                        end
                        
                        shj = [1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,6,6,6,7,7,7];
                        tww = 3;
                        KGEwxy = [];
                        no3sum = [];
                        ftot = [];
                        for p = 0
                            ffnmx=['/taiga/ML/AI_20/Basin' num2str(gid) 'outputs/' rainevent '_SQ/' rainevent '/streamflowarea'];
                            fid=fopen(ffnmx,'rb','ieee-le');
                            tmpdata=fread(fid,inf,'single');
                            fclose(fid);
                            NODAc2o5=reshape(tmpdata,c,r,720);
                            flownoda5 = [];
                            for k = 1:720
                                flownoda5(k) = abs(NODAc2o5(cg,rg,k));
                            end
                            
                            no3 = flownoda5;
                            
                            
                            ffnm=['/taiga/ML/AI_20/Basin' num2str(gid) 'outputs/' rainevent '_SQ/' rainevent '/precipitation'];
                            fid=fopen(ffnm,'rb','ieee-le');
                            tmpdata=fread(fid,inf,'single');
                            fclose(fid);
                            NODAc2o5=reshape(tmpdata,c,r,720);
                            flownoda5 = [];
                            meanr250 = 0;
                            tmpsum = zeros(c,r);
                            for k = 1:720
                                tmp = NODAc2o5(:,:,k);
                                tmpsum = tmpsum+tmp;
                                flownoda5(k) = mean(tmp(intt));
                                meanr250 = meanr250+flownoda5(k);
                            end
                            no2 = flownoda5;%GREEN
                            rainsd = std(tmpsum(intt));
                            
                            cal_window = (stat_left_bd-1)*8+1:(upcap-1)*8+1;
                            
                            if p==0
                                tstmp = find(flownoda5>0.1);% only consider rain>0.1mm/h in the original field
                                % ircicc rain>0.1 hours are a little bit
                                % more than original rain>0.1 hours have to
                                % be consistent for the number of hours for
                                % error modeling
                                rts1 = intersect(tstmp,cal_window);
                                rain00_sumandes{bid,eid} = NODAc2o5(:,:,rts1);
                                eval_indandes{bid,eid} = rts1;
                                simu_floodandes{bid,eid} = no3(rts1);
                                if isempty(rts1)
                                    cvfloodandes(bid,eid) = nan;
                                    maxfloodandes(bid,eid) = nan;
                                    meanfloodandes(bid,eid) = nan;
                                else
                                    cvfloodandes(bid,eid) = mean(no3(rts1))/std(no3(rts1));
                                    maxfloodandes(bid,eid) = max(no3(rts1));
                                    meanfloodandes(bid,eid) = mean(no3(rts1));
                                end
                            end
                            
                            % obviously IRCICC is not done for these
                            % events, then below
                            rain34_sumandes{bid,eid} = nan(c,r,length(rts1));
                            
                            
                        end
                        
                    end
                end
            end
            [bid,eid]
        end
        
    end
    bid
end



%%
% dssumandes = {};
% dosumandes = {};
% pmeansandes = {};
% pchangesandes = {};
% monthnumandes = [];
% yearnumandes = [];

for andis = 1:length(adsall)
    bid = adsall(andis);
    % if bid>404
    %     continue
    % end
    gid = WORLD_events(bid,1);
    tmp = load(['/shared/dondo/home/ml423/DA/Critical matrices/World_Basins/inttind/Basin' num2str(gid) '_ind.mat'],'ind');
    ind = tmp.ind;
    tmp = load(['/shared/dondo/home/ml423/DA/Critical matrices/World_Basins/inttind/Basin' num2str(gid) '_intt.mat'],'intt');
    intt = tmp.intt;       

    d_s = load(['/shared/dondo/home/ml423/DA/Critical matrices/World_Basins/distance/Basin' num2str(gid) 'dis_stream.mat']);
    d_s = d_s.dis_stream;
    d_o = load(['/shared/dondo/home/ml423/DA/Critical matrices/World_Basins/distance/Basin' num2str(gid) 'dis_outlet.mat']);
    d_o = d_o.dis_outlet;
    d_s(ind) = nan;d_o(ind) = nan;
    d_s = d_s+1;%dis_stream is set to be 1 for minimum
    d_s = d_s/max(d_s(:));
    d_o = d_o/max(d_o(:));

    for eid = 1:149
        if isempty(rain00_sumandes{bid,eid})||isempty(rain34_sumandes{bid,eid})
            continue     
        end
        tmp = rain00_sumandes{bid,eid};
        tmp34 = rain34_sumandes{bid,eid};
        tmpdsum = [];
        tmpdoum = [];
        pmeansum = [];
        pchangesum = [];

        for i = 1:size(tmp,3)
           tmps = tmp(:,:,i);
           tmps(ind) = nan;
           numeri = tmps.*d_s;
           tmpds = (sum(numeri(intt)))/sum(tmps(intt));
           tmpdsum = [tmpdsum,tmpds];
           
           numeri = tmps.*d_o;
           tmpdo = (sum(numeri(intt)))/sum(tmps(intt));
           tmpdoum = [tmpdoum,tmpdo];
           
           pmean = mean(tmps(intt));
           pmeansum = [pmeansum,pmean];
           
           tmps34 = tmp34(:,:,i);
           tmpdif = tmps34-tmps;
           
           pchange = mean(tmpdif(intt));
           pchangesum = [pchangesum,pchange];

        end
        dssumandes{bid,eid} = tmpdsum;
        dosumandes{bid,eid} = tmpdoum;
        pmeansandes{bid,eid} = pmeansum;
        pchangesandes{bid,eid} = pchangesum;
        monthnumandes(bid,eid) = mod(floor(WORLD_events(bid,eid+1)/100),100);
        yearnumandes(bid,eid) = floor(WORLD_events(bid,eid+1)/10000);
    end
end


dsallandes = [];
doallandes = [];
pallandes = [];
demallandes = [];
sallandes = [];mallandes = [];yallandes = [];
aallandes = [];
cvfloodallandes = [];
floodmeanallandes = [];floodmaxallandes = [];
floodallandes = [];
targetallandes = [];

toodryproblem = [171:176 193:204];%189-192 ok 20250613

for andis = 1:length(adsall)
    bid = adsall(andis);



    for eid = 1:149
        if ~isempty(dssumandes{bid,eid})&&cvfloodandes(bid,eid)<1000&&~ismember(bid,[toodryproblem])
            numhours = length(dssumandes{bid,eid});
            
            dsallandes = [dsallandes;dssumandes{bid,eid}'];
            doallandes = [doallandes;dosumandes{bid,eid}'];
            pallandes = [pallandes;pmeansandes{bid,eid}'];
            sallandes = [sallandes;slopemeansumandes(bid)*ones(numhours,1)];
            demallandes = [demallandes;demsumandes(bid)*ones(numhours,1)];
            mallandes = [mallandes;monthnumandes(bid,eid)*ones(numhours,1)];
            yallandes = [yallandes;yearnumandes(bid,eid)*ones(numhours,1)];
            aallandes = [aallandes;d_areaandes(bid)*ones(numhours,1)];
            cvfloodallandes = [cvfloodallandes;cvfloodandes(bid,eid)*ones(numhours,1)];
            floodmeanallandes = [floodmeanallandes;meanfloodandes(bid,eid)*ones(numhours,1)];
            floodmaxallandes = [floodmaxallandes;maxfloodandes(bid,eid)*ones(numhours,1)];
            floodallandes = [floodallandes;simu_floodandes{bid,eid}'];
            targetallandes = [targetallandes;pchangesandes{bid,eid}'];
            
        end
        
    end
end
allinandes = [dsallandes';doallandes';pallandes';sallandes';aallandes';demallandes';cvfloodallandes';floodmaxallandes';floodmeanallandes';floodallandes'];
infnum = find(allinandes(7,:)>1000);% a couple of events have 9999 streamflow, therefore rainfall is not reliable for these.
allinandes(:,infnum) = [];
targetallandes(infnum) = [];

%%
% extremedeltap = find(targetallandes>10|targetallandes<-10);
% allinandes(:,extremedeltap) = [];
% targetallandes(extremedeltap) = [];
predic1andes = net1(allinandes);% example case
predic2andes = net2(allinandes);% example case
predic3andes = net3(allinandes);% example case
predic4andes = net4(allinandes);% example case
predic5andes = net5(allinandes);% example case
predic6andes = net6(allinandes);% example case
predic7andes = net7(allinandes);% example case
predic8andes = net8(allinandes);% example case
predic9andes = net9(allinandes);% example case
predic10andes = net10(allinandes);% example case
predic11andes = net11(allinandes);% example case
predic12andes = net12(allinandes);% example case
predic13andes = net13(allinandes);% example case
predic14andes = net14(allinandes);% example case
predic15andes = net15(allinandes);% example case
predic16andes = net16(allinandes);% example case
predic17andes = net17(allinandes);% example case
predic18andes = net18(allinandes);% example case
predic19andes = net19(allinandes);% example case
predic20andes = net20(allinandes);% example case
figure
scatter(predic3,predic11)
figure
scatter(predic3andes,predic4andes)
figure
scatter(targetallandes,predic4andes)
title('predic3')
p1error = find(predic1andes>50);
p3error = find(predic3andes<-13);

p6error = find(predic6andes>100);
p8error = find(predic8andes>100);
loat2 = allinandes(:,p3error)

figure
histogram(targetallandes,[-30:2:30])
xlim([-20 20])
%% 
% apply NN to andes
accuhour = 0;
maxrep = 20;
for andis = 1:length(adsall)
    bid = adsall(andis);
    gid = WORLD_events(bid,1);
    gloc = find(Basins_sum(:,1)==gid);
    
    r = Basins_sum(gloc,6);
    c = Basins_sum(gloc,7);

    if ~ismember(bid,problem_basins)&&~ismember(bid,toodryproblem)
        for eid = 1:149
            if ~isnan(WORLD_events(bid,eid+1))
                event = WORLD_events(bid,eid+1);
                
                ffnm=['/shared/dondo/home/ml423/world_mts/Precip/Event/Basin' num2str(gid) '_input_ERA5_900m_FIX/' num2str(event) '/Precipitation_ERA5_land_oriunit'];
                if exist(ffnm)
                    orir = ['/shared/dondo/home/ml423/world_mts/Precip/Event/Basin' num2str(gid) '_input_ERA5_900m_FIX/' num2str(event) '/'];
                    out_dir2 = ['/shared/dondo/home/ml423/world_mts/Basin' num2str(gid) '_input_ERA5_900m_tmp/' num2str(event) '/'];
                   
                    copyfile([orir '*'],out_dir2)
                    disp('transfer done')
                else
                    disp('error')
                    
                end
                
                if ~isempty(dssumandes{bid,eid})&&cvfloodandes(bid,eid)<1000
                    ffnm=['/shared/dondo/home/ml423/world_mts/Basin' num2str(gid) '_input_ERA5_900m_tmp/' num2str(event) '/Precipitation_ERA5_land_oriunit'];
                    fid=fopen(ffnm,'rb','ieee-le');
                    tmpdata=fread(fid,inf,'single');
                    fclose(fid);
                    NODAc2o5=reshape(tmpdata,c,r,720);
                    
                    kind = eval_indandes{bid,eid};

                    for replic = 1:maxrep
                        predics = eval(['predic' num2str(replic) 'andes']);
             
                        reviserain = [];

                        
                        for k = 1:720
                            tmpori = NODAc2o5(:,:,k);
                            if ismember(k,kind)
                                accuhour = accuhour+1;
                                if replic>1&&k==kind(1)
                                    accuhour = accuhour-length(kind);
                                end
                                tmpupdate = tmpori+predics(accuhour);
                                reviserain(:,:,k) = tmpupdate;
                            else
                                reviserain(:,:,k) = tmpori;
                            end
                            %[k,replic,accuhour]
                        end
                        %bid
                        %eid
                        accuhour
                        bid
                        out_dirp = ['/shared/dondo/home/ml423/world_mts/Basin' num2str(gid) '_input_ERA5_900m_tmp/' num2str(event) '/'];
                        mkdir(out_dirp)
                        
                        reviserain(reviserain<0) = 0;
                        fnm = [out_dirp 'Precipitation_ERA5_land_ori_' num2str(replic,'%3.3d')];
                        raindata = reshape(reviserain,[c,r,720]);
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
                end
                
            end
        end
        
    end
end
%%
% apply Random Forest Regressor to Andes basins
accuhour = 0;
maxrep = 20;
andesRFIn = allinandes';
andesRFIn = andesRFIn(:,[1:5 7 10]);

predic1andes = predict(Mdl_reg1, andesRFIn);
predic2andes = predict(Mdl_reg2, andesRFIn);
predic3andes = predict(Mdl_reg3, andesRFIn);
predic4andes = predict(Mdl_reg4, andesRFIn);
predic5andes = predict(Mdl_reg5, andesRFIn);
predic6andes = predict(Mdl_reg6, andesRFIn);
predic7andes = predict(Mdl_reg7, andesRFIn);
predic8andes = predict(Mdl_reg8, andesRFIn);
predic9andes = predict(Mdl_reg9, andesRFIn);
predic10andes = predict(Mdl_reg10, andesRFIn);
predic11andes = predict(Mdl_reg11, andesRFIn);
predic12andes = predict(Mdl_reg12, andesRFIn);
predic13andes = predict(Mdl_reg13, andesRFIn);
predic14andes = predict(Mdl_reg14, andesRFIn);
predic15andes = predict(Mdl_reg15, andesRFIn);
predic16andes = predict(Mdl_reg16, andesRFIn);
predic17andes = predict(Mdl_reg17, andesRFIn);
predic18andes = predict(Mdl_reg18, andesRFIn);
predic19andes = predict(Mdl_reg19, andesRFIn);
predic20andes = predict(Mdl_reg20, andesRFIn);
for andis = 1:length(adsall)
    bid = adsall(andis);
    gid = WORLD_events(bid,1);
    gloc = find(Basins_sum(:,1)==gid);
    
    r = Basins_sum(gloc,6);
    c = Basins_sum(gloc,7);

    if ~ismember(bid,problem_basins)&&~ismember(bid,toodryproblem)
        for eid = 1:149
            if ~isnan(WORLD_events(bid,eid+1))
                event = WORLD_events(bid,eid+1);
                
                ffnm=['/shared/dondo/home/ml423/world_mts/Precip/Event/Basin' num2str(gid) '_input_ERA5_900m_FIX/' num2str(event) '/Precipitation_ERA5_land_oriunit'];
                if exist(ffnm)
                    orir = ['/shared/dondo/home/ml423/world_mts/Precip/Event/Basin' num2str(gid) '_input_ERA5_900m_FIX/' num2str(event) '/'];
                    out_dir2 = ['/shared/dondo/home/ml423/world_mts/Basin' num2str(gid) '_input_ERA5_900m_tmp/' num2str(event) '/'];
                   
                    copyfile([orir '*'],out_dir2)
                    disp('transfer done')
                else
                    disp('error')
                    
                end
                
                if ~isempty(dssumandes{bid,eid})&&cvfloodandes(bid,eid)<1000
                    ffnm=['/shared/dondo/home/ml423/world_mts/Basin' num2str(gid) '_input_ERA5_900m_tmp/' num2str(event) '/Precipitation_ERA5_land_oriunit'];
                    fid=fopen(ffnm,'rb','ieee-le');
                    tmpdata=fread(fid,inf,'single');
                    fclose(fid);
                    NODAc2o5=reshape(tmpdata,c,r,720);
                    
                    kind = eval_indandes{bid,eid};

                    for replic = 1:maxrep
                        predics = eval(['predic' num2str(replic) 'andes']);
             
                        reviserain = [];

                        
                        for k = 1:720
                            tmpori = NODAc2o5(:,:,k);
                            if ismember(k,kind)
                                accuhour = accuhour+1;
                                if replic>1&&k==kind(1)
                                    accuhour = accuhour-length(kind);
                                end
                                tmpupdate = tmpori+predics(accuhour);
                                reviserain(:,:,k) = tmpupdate;
                            else
                                reviserain(:,:,k) = tmpori;
                            end
                            %[k,replic,accuhour]
                        end
                        %bid
                        %eid
                        accuhour
                        bid
                        out_dirp = ['/shared/dondo/home/ml423/world_mts/Basin' num2str(gid) '_input_ERA5_900m_tmp/' num2str(event) '/'];
                        mkdir(out_dirp)
                        
                        reviserain(reviserain<0) = 0;
                        fnm = [out_dirp 'Precipitation_ERA5_lRFR_ori_' num2str(replic,'%3.3d')];
                        raindata = reshape(reviserain,[c,r,720]);
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
                end
                
            end
        end
        
    end
end

%% Apply Alps model to basins in Aisa
% calculate slopes
himaall = find(gauge_events_first_last(:,5)<1000);
hiall = [];
for ki = 1:length(himaall)
    bid = find(WORLD_events(:,1)==gauge_events_first_last(himaall(ki),1));
    if bid>411
            hiall = [hiall;bid];
    end
end
for hima = 1:length(hiall)
    bid = hiall(hima);
    gid = WORLD_events(bid,1);
    gloc = find(Basins_sum(:,1)==gid);
    
    r = Basins_sum(gloc,6);
    c = Basins_sum(gloc,7);
    tmp = load(['/shared/dondo/home/ml423/DA/Critical matrices/World_Basins/inttind/Basin' num2str(gid) '_ind.mat'],'ind');
    ind = tmp.ind;
    tmp = load(['/shared/dondo/home/ml423/DA/Critical matrices/World_Basins/inttind/Basin' num2str(gid) '_intt.mat'],'intt');
    intt = tmp.intt; 
    ffnm=['/shared/dondo/home/ml423/world_mts/Basin' num2str(gid) '_fixed_900m/dem.bin'];
    fid=fopen(ffnm,'rb','ieee-le');
    tmpdata=fread(fid,inf,'single');
    fclose(fid);
    dem = reshape(tmpdata,c,r);
    %dem(ind) = nan;
    demsumhima(bid,1) = mean(dem(intt));
    
    ffnm=['/shared/dondo/home/ml423/world_mts/Basin' num2str(gid) '_fixed_900m/fdr.bin'];
    fid=fopen(ffnm);
    tmpdata=fread(fid,'ubit16');
    fclose(fid);
    fdr = reshape(tmpdata,c,r);
    slopepix = nan(c,r);
%     if bid == 341
%         excep=3;
%     else
%         excep=2;
%     end
    for i = excep:c-1
        for j = 2:r-1
%             if ismember(c*j+i,ind)
%                break 
%             end
            tmpz = fdr(i,j);
            tmp1 = i;tmp2 = j;
            d1 = dem(i,j);

            if tmpz == 1
                ik = tmp1+1;
                jk = tmp2-1;
                dis = 900*sqrt(2);
            elseif  tmpz == 2
                ik = tmp1+1;
                jk = tmp2;
                dis = 900;
            elseif  tmpz == 4
                ik = tmp1+1;
                jk = tmp2+1;
                dis = 900*sqrt(2);
            elseif  tmpz == 8
                ik = tmp1;
                jk = tmp2+1;
                dis = 900;
            elseif  tmpz == 16
                ik = tmp1-1;
                jk = tmp2+1;
                dis = 900*sqrt(2);
            elseif  tmpz == 32
                ik = tmp1-1;
                jk = tmp2;
                dis = 900;
            elseif  tmpz == 64
                ik = tmp1-1;
                jk = tmp2-1;
                dis = 900*sqrt(2);
            elseif  tmpz == 128
                ik = tmp1;
                jk = tmp2-1;
                dis = 900;
            end
            d2 = dem(ik,jk);
            slopepix(i,j) = (d1-d2)/dis;
            if slopepix(i,j)==0
                slopepix(i,j)=0.0001;
            end

        end
    end
    slopemeansumhima(bid,1) = mean(slopepix(intt));
end

for hima = 1:length(hiall)
    bid = hiall(hima);
    gid = WORLD_events(bid,1);
    
    gloc = find(Basins_sum(:,1)==gid);
    
    r = Basins_sum(gloc,6);
    c = Basins_sum(gloc,7);
    rg = Basins_sum(gloc,8);
    cg = Basins_sum(gloc,9);
    d_areahima(bid) = Basins_sum(gloc,5);
    
    tmp = load(['/shared/dondo/home/ml423/DA/Critical matrices/World_Basins/inttind/Basin' num2str(gid) '_ind.mat'],'ind');
    ind = tmp.ind;
    tmp = load(['/shared/dondo/home/ml423/DA/Critical matrices/World_Basins/inttind/Basin' num2str(gid) '_intt.mat'],'intt');
    intt = tmp.intt;
    ffnm=['/shared/dondo/home/ml423/world_mts/Basin' num2str(gid) '_fixed_900m/dem.bin'];
    fid=fopen(ffnm,'rb','ieee-le');
    tmpdata=fread(fid,inf,'single');
    fclose(fid);
    dem = reshape(tmpdata,c,r);
    dem(ind) = nan;
    
    ixx = find(~isnan(WORLD_events(bid,:)));
    etot = length(ixx)-1;
    if etot == 0
        continue
    else
        
        for eid = 1:etot
            close all
            if isnan(event_quality(bid,eid+1))
                disp('Event quality is bad, peak too early or too late')
                continue
            else
                rainevent = num2str(WORLD_events(bid,eid+1));
                yth = floor(WORLD_events(bid,eid+1)/10000);
                dth = mod(WORLD_events(bid,eid+1),100);
                mth = floor((WORLD_events(bid,eid+1)-yth*10000)/100);
                
                thresdate = datetime(yth,mth,dth);
                
                durati = (thresdate-days(15)):days(1):(thresdate+days(15));
                [dury,durm,durd] = ymd(durati);
                tst1 = num2str(durm(1)*100+durd(1),'%4.4d');
                tst2 = num2str(durm(30)*100+durd(30),'%4.4d');
                
                % the condition below has to change, most of events in
                % hima is not used IRCICC because no computational resouce
                ffnm=['/taiga/ML/AI_20/Basin' num2str(gid) 'outputs/' rainevent '_SQ/' rainevent '/streamflowarea'];
                fid=fopen(ffnm,'rb','ieee-le');
                if fid == -1
                    continue
                else
                    fclose(fid);
                    
                    ffnm = ['/shared/dondo/home/ml423/world_mts_runs/Basins_events_chars/'];
                    fid = fopen([ffnm 'Basin' num2str(gid) '_event' rainevent '.txt'],'r');
                    if fid == -1
                        fclose(fid);
                        continue
                    end
                    windowtimes = fscanf(fid,'%d\n',[4,1]);
                    fclose(fid);
                    if length(windowtimes)<4
                        continue
                    end
                    %close all
                    
                    
                    if 1==1
                        
                        %fh = APL_frontcut(bid,eid);
                        
                        tstr = [num2str(dury(1)) '-' num2str(durm(1),'%2.2d') '/' num2str(durd(1),'%2.2d') '-' num2str(durm(30),'%2.2d') '/' num2str(durd(30),'%2.2d')];

                        obstep = 90;
                        
                        tmp = load(['/shared/dondo/home/ml423/world_mts/obs_events/events/obs_' num2str(gid) '_event' rainevent '.mat']);
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
                        
                        
                        
                        tmp = load(['/shared/dondo/home/ml423/world_mts/Precip/Event/ERA5_Basin' num2str(gid) '_' rainevent '.mat']);
                        s4dbkc = tmp.rain;
                        
                        MCp = zeros(c,r);
                        pts = [];
                        dbkcr = 0;rmax = [];
                        for i = 1:720
                            MCp = MCp + s4dbkc{i};
                            pts(i) = mean(s4dbkc{i}(intt));
                            dbkcr = dbkcr+pts(i);
                            rmax(i) = max(s4dbkc{i}(intt));
                        end
                        MCp(ind) = nan;
                        ori_mean = nanmean(MCp(:));
                        ama = jet(40);ama(1,:) = [1 1 1];
                        
                        
                        % for Nature Science paper, consider the major events only
                        
                        launchp = find(rmax>0.5);% when rainfall >1mm/hr launch particles
                        ptmaxthre = 0.5*max(rmax);% used as a threshold to separate events
                        la2 = find(launchp<(windowtimes(4)*8));
                        la1 = find(launchp<(windowtimes(3)*8));
                        accu0 = 0;
                        if isempty(la1)
                            continue
                        end
                        for tmpla = launchp(la1(end)):-1:1
                            if ismember(tmpla,launchp)
                                accu0 = 0;
                                rainbegintime = tmpla;
                            else
                                accu0 = accu0+1;
                            end
                            
                            if area>1000
                                noraininterval = 12;% for big basins, it is likely there is always rain in the basin, therefore using a small norain limit to distinguish different rain events.
                            else
                                noraininterval = 24;
                            end
                            if accu0>noraininterval
                                % consecutive x hours of basin max rain less than 0.5mm/hr,
                                % this should be counted as a different rain cell, should not
                                % be tracked!
                                % one exception is if another major rain is 48-24 ahead, should be
                                % considered as well
                                rsmallest = max([(tmpla-noraininterval),1]);
                                
                                reconsider_period = rsmallest:tmpla;
                                %reconsider_period = (tmpla-noraininterval):tmpla;
                                secondary_rain = find(rmax(reconsider_period)>0.5);
                                
                                if max(rmax(reconsider_period))>ptmaxthre
                                    % adjust the rainbegintime here
                                    rainbegintime = reconsider_period(secondary_rain(1));
                                end
                                smalltmpla = max([floor(tmpla/8),2]);
                                if obsnow(smalltmpla)>obsnow(smalltmpla-1)
                                    % meaning obsflow is still rising.
                                    if ~isempty(secondary_rain)
                                        rainbegintime = reconsider_period(secondary_rain(1));
                                    end
                                end
                                
                                break
                            end
                        end
                        
                        la0 = find(launchp == rainbegintime);
                        
                        refinedlaunchp = launchp(la0:la2(end));
                        stat_left_bd_hydro = windowtimes(1);% systematic defined
                        stat_left_bd_rain = floor(rainbegintime/8);%target rainevent adjusted
                        stat_left_bd = max([stat_left_bd_hydro,stat_left_bd_rain]);
                        
                        upcap = min([90 (windowtimes(4)+6)]);
                        
                        if upcap<=stat_left_bd
                            continue
                        end
                        
                        shj = [1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,6,6,6,7,7,7];
                        tww = 3;
                        KGEwxy = [];
                        no3sum = [];
                        ftot = [];
                        for p = 0
                            ffnmx=['/taiga/ML/AI_20/Basin' num2str(gid) 'outputs/' rainevent '_SQ/' rainevent '/streamflowarea'];
                            fid=fopen(ffnmx,'rb','ieee-le');
                            tmpdata=fread(fid,inf,'single');
                            fclose(fid);
                            NODAc2o5=reshape(tmpdata,c,r,720);
                            flownoda5 = [];
                            for k = 1:720
                                flownoda5(k) = abs(NODAc2o5(cg,rg,k));
                            end
                            
                            no3 = flownoda5;
                            
                            
                            ffnm=['/taiga/ML/AI_20/Basin' num2str(gid) 'outputs/' rainevent '_SQ/' rainevent '/precipitation'];
                            fid=fopen(ffnm,'rb','ieee-le');
                            tmpdata=fread(fid,inf,'single');
                            fclose(fid);
                            NODAc2o5=reshape(tmpdata,c,r,720);
                            flownoda5 = [];
                            meanr250 = 0;
                            tmpsum = zeros(c,r);
                            for k = 1:720
                                tmp = NODAc2o5(:,:,k);
                                tmpsum = tmpsum+tmp;
                                flownoda5(k) = mean(tmp(intt));
                                meanr250 = meanr250+flownoda5(k);
                            end
                            no2 = flownoda5;%GREEN
                            rainsd = std(tmpsum(intt));
                            
                            cal_window = (stat_left_bd-1)*8+1:(upcap-1)*8+1;
                            
                            if p==0
                                tstmp = find(flownoda5>0.1);% only consider rain>0.1mm/h in the original field
                                % ircicc rain>0.1 hours are a little bit
                                % more than original rain>0.1 hours have to
                                % be consistent for the number of hours for
                                % error modeling
                                rts1 = intersect(tstmp,cal_window);
                                rain00_sumhima{bid,eid} = NODAc2o5(:,:,rts1);
                                eval_indhima{bid,eid} = rts1;
                                simu_floodhima{bid,eid} = no3(rts1);
                                if isempty(rts1)
                                    cvfloodhima(bid,eid) = nan;
                                    maxfloodhima(bid,eid) = nan;
                                    meanfloodhima(bid,eid) = nan;
                                else
                                    cvfloodhima(bid,eid) = mean(no3(rts1))/std(no3(rts1));
                                    maxfloodhima(bid,eid) = max(no3(rts1));
                                    meanfloodhima(bid,eid) = mean(no3(rts1));
                                end
                            end
                            
                            rain34_sumhima{bid,eid} = nan(c,r,length(rts1));
                            
                            
                        end
                        
                    end
                end
            end
            [bid,eid]
        end
        
    end
    bid
end

%
for hima = 1:length(hiall)
    bid = hiall(hima);
    gid = WORLD_events(bid,1);
    tmp = load(['/shared/dondo/home/ml423/DA/Critical matrices/World_Basins/inttind/Basin' num2str(gid) '_ind.mat'],'ind');
    ind = tmp.ind;
    tmp = load(['/shared/dondo/home/ml423/DA/Critical matrices/World_Basins/inttind/Basin' num2str(gid) '_intt.mat'],'intt');
    intt = tmp.intt;       

    d_s = load(['/shared/dondo/home/ml423/DA/Critical matrices/World_Basins/distance/Basin' num2str(gid) 'dis_stream.mat']);
    d_s = d_s.dis_stream;
    d_o = load(['/shared/dondo/home/ml423/DA/Critical matrices/World_Basins/distance/Basin' num2str(gid) 'dis_outlet.mat']);
    d_o = d_o.dis_outlet;
    d_s(ind) = nan;d_o(ind) = nan;
    d_s = d_s+1;%dis_stream is set to be 1 for minimum
    d_s = d_s/max(d_s(:));
    d_o = d_o/max(d_o(:));

    for eid = 1:58
        if isempty(rain00_sumhima{bid,eid})||isempty(rain34_sumhima{bid,eid})
            continue     
        end
        tmp = rain00_sumhima{bid,eid};
        tmp34 = rain34_sumhima{bid,eid};
        tmpdsum = [];
        tmpdoum = [];
        pmeansum = [];
        pchangesum = [];

        for i = 1:size(tmp,3)
           tmps = tmp(:,:,i);
           tmps(ind) = nan;
           numeri = tmps.*d_s;
           tmpds = (sum(numeri(intt)))/sum(tmps(intt));
           tmpdsum = [tmpdsum,tmpds];
           
           numeri = tmps.*d_o;
           tmpdo = (sum(numeri(intt)))/sum(tmps(intt));
           tmpdoum = [tmpdoum,tmpdo];
           
           pmean = mean(tmps(intt));
           pmeansum = [pmeansum,pmean];
           
           tmps34 = tmp34(:,:,i);
           tmpdif = tmps34-tmps;
           
           pchange = mean(tmpdif(intt));
           pchangesum = [pchangesum,pchange];

        end
        dssumhima{bid,eid} = tmpdsum;
        dosumhima{bid,eid} = tmpdoum;
        pmeanshima{bid,eid} = pmeansum;
        pchangeshima{bid,eid} = pchangesum;
        monthnumhima(bid,eid) = mod(floor(WORLD_events(bid,eid+1)/100),100);
        yearnumhima(bid,eid) = floor(WORLD_events(bid,eid+1)/10000);
    end
end


dsallhima = [];
doallhima = [];
pallhima = [];
demallhima = [];
sallhima = [];mallhima = [];yallhima = [];
aallhima = [];
cvfloodallhima = [];
floodmeanallhima = [];floodmaxallhima = [];
floodallhima = [];
targetallhima = [];

toodryproblem = [];

for hima = 1:length(hiall)
    bid = hiall(hima);
    gid = WORLD_events(bid,1);


    for eid = 1:58
        if ~isempty(dssumhima{bid,eid})&&cvfloodhima(bid,eid)<1000&&~ismember(bid,[toodryproblem])
            numhours = length(dssumhima{bid,eid});
            
            dsallhima = [dsallhima;dssumhima{bid,eid}'];
            doallhima = [doallhima;dosumhima{bid,eid}'];
            pallhima = [pallhima;pmeanshima{bid,eid}'];
            sallhima = [sallhima;slopemeansumhima(bid)*ones(numhours,1)];
            demallhima = [demallhima;demsumhima(bid)*ones(numhours,1)];
            mallhima = [mallhima;monthnumhima(bid,eid)*ones(numhours,1)];
            yallhima = [yallhima;yearnumhima(bid,eid)*ones(numhours,1)];
            aallhima = [aallhima;d_areahima(bid)*ones(numhours,1)];
            cvfloodallhima = [cvfloodallhima;cvfloodhima(bid,eid)*ones(numhours,1)];
            floodmeanallhima = [floodmeanallhima;meanfloodhima(bid,eid)*ones(numhours,1)];
            floodmaxallhima = [floodmaxallhima;maxfloodhima(bid,eid)*ones(numhours,1)];
            floodallhima = [floodallhima;simu_floodhima{bid,eid}'];
            targetallhima = [targetallhima;pchangeshima{bid,eid}'];
            
        end
        
    end
end

%
allinhima = [dsallhima';doallhima';pallhima';sallhima';aallhima';demallhima';cvfloodallhima';floodmaxallhima';floodmeanallhima';floodallhima'];
infnum = find(allinhima(7,:)>100);% a couple of events have 9999 streamflow, therefore rainfall is not reliable for these.
allinhima(:,infnum) = [];
targetallhima(infnum) = [];

% extremedeltap = find(targetallhima>10|targetallhima<-10);
% allinhima(:,extremedeltap) = [];
% targetallhima(extremedeltap) = [];
predic1hima = net1(allinhima);% example case
predic2hima = net2(allinhima);% example case
predic3hima = net3(allinhima);% example case
predic4hima = net4(allinhima);% example case
predic5hima = net5(allinhima);% example case
predic6hima = net6(allinhima);% example case
predic7hima = net7(allinhima);% example case
predic8hima = net8(allinhima);% example case
predic9hima = net9(allinhima);% example case
predic10hima = net10(allinhima);% example case
predic11hima = net11(allinhima);% example case
predic12hima = net12(allinhima);% example case
predic13hima = net13(allinhima);% example case
predic14hima = net14(allinhima);% example case
predic15hima = net15(allinhima);% example case
predic16hima = net16(allinhima);% example case
predic17hima = net17(allinhima);% example case
predic18hima = net18(allinhima);% example case
predic19hima = net19(allinhima);% example case
predic20hima = net20(allinhima);% example case
% figure
% scatter(predic3,predic11)
% figure
% scatter(predic3hima,predic4hima)
% figure
% scatter(targetallhima,predic4hima)
% title('predic3')
% p1error = find(predic1hima>50);
% p3error = find(predic3hima<-13);
% 
% p6error = find(predic6hima>100);
% p8error = find(predic8hima>100);
% loat2 = allinhima(:,p3error)
% 
% figure
% histogram(targetallhima,[-30:2:30])
% xlim([-20 20])
%
accuhour = 0;
for hima = 1:length(hiall)
    bid = hiall(hima);
    gid = WORLD_events(bid,1);
    gloc = find(Basins_sum(:,1)==gid);
    
    r = Basins_sum(gloc,6);
    c = Basins_sum(gloc,7);

    if ~ismember(bid,problem_basins)&&~ismember(bid,toodryproblem)
        for eid = 1:149
            if ~isnan(WORLD_events(bid,eid+1))
                event = WORLD_events(bid,eid+1);
                
                ffnm=['/shared/dondo/home/ml423/world_mts/Precip/Event/Basin' num2str(gid) '_input_ERA5_900m_FIX/' num2str(event) '/Precipitation_ERA5_land_oriunit'];
                if exist(ffnm)
                    orir = ['/shared/dondo/home/ml423/world_mts/Precip/Event/Basin' num2str(gid) '_input_ERA5_900m_FIX/' num2str(event) '/'];
                    out_dir2 = ['/shared/dondo/home/ml423/world_mts/Basin' num2str(gid) '_input_ERA5_900m_tmp/' num2str(event) '/'];
                   
                    copyfile([orir '*'],out_dir2)
                    disp('transfer done')
                else
                    disp('error')
                    
                end
                
                if ~isempty(dssumhima{bid,eid})&&cvfloodhima(bid,eid)<1000
                    ffnm=['/shared/dondo/home/ml423/world_mts/Basin' num2str(gid) '_input_ERA5_900m_tmp/' num2str(event) '/Precipitation_ERA5_land_oriunit'];
                    fid=fopen(ffnm,'rb','ieee-le');
                    tmpdata=fread(fid,inf,'single');
                    fclose(fid);
                    NODAc2o5=reshape(tmpdata,c,r,720);
                    
                    kind = eval_indhima{bid,eid};

                    for replic = 1:20
                        predics = eval(['predic' num2str(replic) 'hima']);
             
                        reviserain = [];

                        
                        for k = 1:720
                            tmpori = NODAc2o5(:,:,k);
                            if ismember(k,kind)
                                accuhour = accuhour+1;
                                if replic>1&&k==kind(1)
                                    accuhour = accuhour-length(kind);
                                end
                                tmpupdate = tmpori+predics(accuhour);
                                reviserain(:,:,k) = tmpupdate;
                            else
                                reviserain(:,:,k) = tmpori;
                            end
                            %[k,replic,accuhour]
                        end
                        %bid
                        %eid
                        accuhour
                        out_dirp = ['/shared/dondo/home/ml423/world_mts/Basin' num2str(gid) '_input_ERA5_900m_tmp/' num2str(event) '/'];
                        mkdir(out_dirp)
                        
                        reviserain(reviserain<0) = 0;
                        fnm = [out_dirp 'Precipitation_ERA5_land_ori_' num2str(replic,'%3.3d')];
                        raindata = reshape(reviserain,[c,r,720]);
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
                end
                
            end
        end
        
    end
end

%% Random forest models, check the impact of initial weights
% conclusion: RF is not sensitive to initial weights. If NN is used, then
% an ensemble models are needed to account for nonlinear optimization
% problems.
X_reg = allin';X_reg = X_reg(:,[1:5 7 10]);
Y_reg = targ';
X_test = allin_test';X_test = X_test(:,[1:5 7 10]);

disp('--- Regression Example ---');
rng('default'); % For reproducibility
rng(45); % For reproducibility, 0 42 43 44 45 
% Train a TreeBagger model for regression
numTrees_reg = 30;
Mdl_reg5 = TreeBagger(numTrees_reg, X_reg, Y_reg, ...
                     'Method', 'regression', ...
                     'OOBPrediction', 'on', ...
                     'OOBPredictorImportance', 'on', ...
                     'MinLeafSize', 10); % Enable out-of-bag error estimation
rng(46); % For reproducibility, 0 42 43 44 45 
numTrees_reg = 30;
Mdl_reg6 = TreeBagger(numTrees_reg, X_reg, Y_reg, ...
                     'Method', 'regression', ...
                     'OOBPrediction', 'on', ...
                     'OOBPredictorImportance', 'on', ...
                     'MinLeafSize', 10); % Enable out-of-bag error estimation
rng(47); % For reproducibility, 0 42 43 44 45 
numTrees_reg = 30;
Mdl_reg7 = TreeBagger(numTrees_reg, X_reg, Y_reg, ...
                     'Method', 'regression', ...
                     'OOBPrediction', 'on', ...
                     'OOBPredictorImportance', 'on', ...
                     'MinLeafSize', 10); % Enable out-of-bag error estimation
rng(48); % For reproducibility, 0 42 43 44 45 
numTrees_reg = 30;
Mdl_reg8 = TreeBagger(numTrees_reg, X_reg, Y_reg, ...
                     'Method', 'regression', ...
                     'OOBPrediction', 'on', ...
                     'OOBPredictorImportance', 'on', ...
                     'MinLeafSize', 10); % Enable out-of-bag error estimation
rng(49); % For reproducibility, 0 42 43 44 45 
numTrees_reg = 30;
Mdl_reg9 = TreeBagger(numTrees_reg, X_reg, Y_reg, ...
                     'Method', 'regression', ...
                     'OOBPrediction', 'on', ...
                     'OOBPredictorImportance', 'on', ...
                     'MinLeafSize', 10); % Enable out-of-bag error estimation
rng(50); % For reproducibility, 0 42 43 44 45 
numTrees_reg = 30;
Mdl_reg10 = TreeBagger(numTrees_reg, X_reg, Y_reg, ...
                     'Method', 'regression', ...
                     'OOBPrediction', 'on', ...
                     'OOBPredictorImportance', 'on', ...
                     'MinLeafSize', 10); % Enable out-of-bag error estimation
rng(51); % For reproducibility, 0 42 43 44 45 
numTrees_reg = 30;
Mdl_reg11 = TreeBagger(numTrees_reg, X_reg, Y_reg, ...
                     'Method', 'regression', ...
                     'OOBPrediction', 'on', ...
                     'OOBPredictorImportance', 'on', ...
                     'MinLeafSize', 10); % Enable out-of-bag error estimation
rng(52); % For reproducibility, 0 42 43 44 45 
numTrees_reg = 30;
Mdl_reg12 = TreeBagger(numTrees_reg, X_reg, Y_reg, ...
                     'Method', 'regression', ...
                     'OOBPrediction', 'on', ...
                     'OOBPredictorImportance', 'on', ...
                     'MinLeafSize', 10); % Enable out-of-bag error estimation
rng(53); % For reproducibility, 0 42 43 44 45 
numTrees_reg = 30;
Mdl_reg13 = TreeBagger(numTrees_reg, X_reg, Y_reg, ...
                     'Method', 'regression', ...
                     'OOBPrediction', 'on', ...
                     'OOBPredictorImportance', 'on', ...
                     'MinLeafSize', 10); % Enable out-of-bag error estimation
rng(54); % For reproducibility, 0 42 43 44 45 
numTrees_reg = 30;
Mdl_reg14 = TreeBagger(numTrees_reg, X_reg, Y_reg, ...
                     'Method', 'regression', ...
                     'OOBPrediction', 'on', ...
                     'OOBPredictorImportance', 'on', ...
                     'MinLeafSize', 10); % Enable out-of-bag error estimation
rng(55); % For reproducibility, 0 42 43 44 45 
numTrees_reg = 30;
Mdl_reg15 = TreeBagger(numTrees_reg, X_reg, Y_reg, ...
                     'Method', 'regression', ...
                     'OOBPrediction', 'on', ...
                     'OOBPredictorImportance', 'on', ...
                     'MinLeafSize', 10); % Enable out-of-bag error estimation
rng(56); % For reproducibility, 0 42 43 44 45 
numTrees_reg = 30;
Mdl_reg16 = TreeBagger(numTrees_reg, X_reg, Y_reg, ...
                     'Method', 'regression', ...
                     'OOBPrediction', 'on', ...
                     'OOBPredictorImportance', 'on', ...
                     'MinLeafSize', 10); % Enable out-of-bag error estimation
rng(57); % For reproducibility, 0 42 43 44 45 
numTrees_reg = 30;
Mdl_reg17 = TreeBagger(numTrees_reg, X_reg, Y_reg, ...
                     'Method', 'regression', ...
                     'OOBPrediction', 'on', ...
                     'OOBPredictorImportance', 'on', ...
                     'MinLeafSize', 10); % Enable out-of-bag error estimation
rng(58); % For reproducibility, 0 42 43 44 45 
numTrees_reg = 30;
Mdl_reg18 = TreeBagger(numTrees_reg, X_reg, Y_reg, ...
                     'Method', 'regression', ...
                     'OOBPrediction', 'on', ...
                     'OOBPredictorImportance', 'on', ...
                     'MinLeafSize', 10); % Enable out-of-bag error estimation
rng(59); % For reproducibility, 0 42 43 44 45 
numTrees_reg = 30;
Mdl_reg19 = TreeBagger(numTrees_reg, X_reg, Y_reg, ...
                     'Method', 'regression', ...
                     'OOBPrediction', 'on', ...
                     'OOBPredictorImportance', 'on', ...
                     'MinLeafSize', 10); % Enable out-of-bag error estimation
rng(60); % For reproducibility, 0 42 43 44 45 
numTrees_reg = 30;
Mdl_reg20 = TreeBagger(numTrees_reg, X_reg, Y_reg, ...
                     'Method', 'regression', ...
                     'OOBPrediction', 'on', ...
                     'OOBPredictorImportance', 'on', ...
                     'MinLeafSize', 10); % Enable out-of-bag error estimation
% save(['/shared/dondo/home/ml423/world_mts_runs/AI_test_basin_event_fold_' num2str(kfold) '_RF' num2str(1,'%2.2d') '.mat'],'Mdl_reg1')
% save(['/shared/dondo/home/ml423/world_mts_runs/AI_test_basin_event_fold_' num2str(kfold) '_RF' num2str(2,'%2.2d') '.mat'],'Mdl_reg2')
% save(['/shared/dondo/home/ml423/world_mts_runs/AI_test_basin_event_fold_' num2str(kfold) '_RF' num2str(3,'%2.2d') '.mat'],'Mdl_reg3')
% save(['/shared/dondo/home/ml423/world_mts_runs/AI_test_basin_event_fold_' num2str(kfold) '_RF' num2str(4,'%2.2d') '.mat'],'Mdl_reg4')
% save(['/shared/dondo/home/ml423/world_mts_runs/AI_test_basin_event_fold_' num2str(kfold) '_RF' num2str(5,'%2.2d') '.mat'],'Mdl_reg5')
% save(['/shared/dondo/home/ml423/world_mts_runs/AI_test_basin_event_fold_' num2str(kfold) '_RF' num2str(6,'%2.2d') '.mat'],'Mdl_reg6')
% save(['/shared/dondo/home/ml423/world_mts_runs/AI_test_basin_event_fold_' num2str(kfold) '_RF' num2str(7,'%2.2d') '.mat'],'Mdl_reg7')
% save(['/shared/dondo/home/ml423/world_mts_runs/AI_test_basin_event_fold_' num2str(kfold) '_RF' num2str(8,'%2.2d') '.mat'],'Mdl_reg8')
% save(['/shared/dondo/home/ml423/world_mts_runs/AI_test_basin_event_fold_' num2str(kfold) '_RF' num2str(9,'%2.2d') '.mat'],'Mdl_reg9')
% save(['/shared/dondo/home/ml423/world_mts_runs/AI_test_basin_event_fold_' num2str(kfold) '_RF' num2str(10,'%2.2d') '.mat'],'Mdl_reg10')
% save(['/shared/dondo/home/ml423/world_mts_runs/AI_test_basin_event_fold_' num2str(kfold) '_RF' num2str(11,'%2.2d') '.mat'],'Mdl_reg11')
% save(['/shared/dondo/home/ml423/world_mts_runs/AI_test_basin_event_fold_' num2str(kfold) '_RF' num2str(12,'%2.2d') '.mat'],'Mdl_reg12')
% save(['/shared/dondo/home/ml423/world_mts_runs/AI_test_basin_event_fold_' num2str(kfold) '_RF' num2str(13,'%2.2d') '.mat'],'Mdl_reg13')
% save(['/shared/dondo/home/ml423/world_mts_runs/AI_test_basin_event_fold_' num2str(kfold) '_RF' num2str(14,'%2.2d') '.mat'],'Mdl_reg14')
% save(['/shared/dondo/home/ml423/world_mts_runs/AI_test_basin_event_fold_' num2str(kfold) '_RF' num2str(15,'%2.2d') '.mat'],'Mdl_reg15')
% save(['/shared/dondo/home/ml423/world_mts_runs/AI_test_basin_event_fold_' num2str(kfold) '_RF' num2str(16,'%2.2d') '.mat'],'Mdl_reg16')
% save(['/shared/dondo/home/ml423/world_mts_runs/AI_test_basin_event_fold_' num2str(kfold) '_RF' num2str(17,'%2.2d') '.mat'],'Mdl_reg17')
% save(['/shared/dondo/home/ml423/world_mts_runs/AI_test_basin_event_fold_' num2str(kfold) '_RF' num2str(18,'%2.2d') '.mat'],'Mdl_reg18')
% save(['/shared/dondo/home/ml423/world_mts_runs/AI_test_basin_event_fold_' num2str(kfold) '_RF' num2str(19,'%2.2d') '.mat'],'Mdl_reg19')
% save(['/shared/dondo/home/ml423/world_mts_runs/AI_test_basin_event_fold_' num2str(kfold) '_RF' num2str(20,'%2.2d') '.mat'],'Mdl_reg20')
% Predict on new data (using the same data for demonstration)
Yfit_test = predict(Mdl_reg5, X_test);

% Evaluate performance (e.g., out-of-bag error)
oobError_reg = oobError(Mdl_reg5);
disp(['Regression OOB Error: ', num2str(oobError_reg(end))]);

% Plotting results for regression (optional)
figure;
scatter(targetall_test, Yfit_test, 'b.');
%hold on;
%plot(Y_test, Yfit_test, 'r.');
xlabel('target');
ylabel('predict');
title('TreeBagger Regression Prediction');
%legend('True MPG', 'Predicted MPG');


%% [key] plot AI outputs for Alps test-f(basin, event) 

clear tmpp
clear tmp
clear predatatmp
clear data3D
addpath('/shared/dondo/home/ml423/HMIOP/NoDAre/New era/')
addpath('/shared/dondo/home/ml423/DA/Critical matrices/')
addpath('/shared/dondo/home/ml423/DA/slanCM/')

tmp = load(['/shared/dondo/home/ml423/DA/Critical matrices/world_selected_gauge_events_allinfo_v3.mat'],'gauge_events_first_last');
Basins_sum = tmp.gauge_events_first_last;clear tmp;
tmp = load(['/shared/dondo/home/ml423/DA/Critical matrices/WORLD_events_v3.mat'],'WORLD_events');
WORLD_events = tmp.WORLD_events;clear tmp;

flowtype = {'streamflowarea','interflowarea','overlandflowarea','baseflowarea'};
ftp = 1;
AItool = 'A15';
alps_AIstats = [];
aiens = 5;
%bes = sortrows(basin_event_test);
bes = sortrows(alleveori);
clear unfinish_alps
%illo = find(bes(:,1)==29);
probAI_alpss = [63 35];
rain_AI = cell(522,149);
for ii = 1:length(bes)
    %bid = [1:129 132:133 466:480 482:504 507:522]
    disp(ii)
    %ii = find(bes(:,1)==probAI_alpss(iii,1)&bes(:,2)==probAI_alpss(iii,2));
    bid = bes(ii,1);

    gid = WORLD_events(bid,1);
    
    gloc = find(Basins_sum(:,1)==gid);
    
    r = Basins_sum(gloc,6);
    c = Basins_sum(gloc,7);
    rg = Basins_sum(gloc,8);
    cg = Basins_sum(gloc,9);
    area = Basins_sum(gloc,5);

    tmp = load(['/shared/dondo/home/ml423/DA/Critical matrices/World_Basins/inttind/Basin' num2str(gid) '_ind.mat'],'ind');
    ind = tmp.ind;
    tmp = load(['/shared/dondo/home/ml423/DA/Critical matrices/World_Basins/inttind/Basin' num2str(gid) '_intt.mat'],'intt');
    intt = tmp.intt;
    %{
    ffnm=['/shared/dondo/home/ml423/world_mts/Basin' num2str(gid) '_fixed_900m/dem.bin'];
    fid=fopen(ffnm,'rb','ieee-le');
    tmpdata=fread(fid,inf,'single');
    fclose(fid);
    dem = reshape(tmpdata,c,r);
    dem(ind) = nan;
    demf = dem';
    colo = slanCM('terrain',80);
    colo(1,:) = [0.8 0.8 0.8];
    figure
    imagesc(demf)
    colormap(colo)
    
    ffnm=['/shared/dondo/home/ml423/world_mts/Basin' num2str(gid) '_fixed_900m/countfile.out'];
    fid=fopen(ffnm);
    facc=fscanf(fid,'%d',[c,r]);
    fclose(fid);
    facc(ind) = nan;
    [streamx,streamy] = find(facc>5);
    %}
    
    
    
    ICn = 4140;
    ICC = 0;% IF ICC is part of the algo
    
    fileind = -100000*ICn;
    KGEsu = [];
    evesu = [];
    budgetalldbkc = [];
    budgetallirc = [];
    
    ixx = find(~isnan(WORLD_events(bid,:)));
    etot = length(ixx)-1;
    if etot == 0
        continue
    else
        if etot>149
           etot = 149 
        end
        for eidd = 1
            eid = bes(ii,2);
            close all
            % the if comments below are data quality control of the Alps
            % and hima
            %if isnan(event_quality(bid,eid+1))||isnan(kge_optm(bid,eid))
            if ~ismember([bid,eid],alleveori,'rows')
                disp('Event quality is bad, peak too early or too late')
                continue
            else
                rainevent = num2str(WORLD_events(bid,eid+1));
                yth = floor(WORLD_events(bid,eid+1)/10000);
                dth = mod(WORLD_events(bid,eid+1),100);
                mth = floor((WORLD_events(bid,eid+1)-yth*10000)/100);

                thresdate = datetime(yth,mth,dth);

                durati = (thresdate-days(15)):days(1):(thresdate+days(15));
                [dury,durm,durd] = ymd(durati);
                tst1 = num2str(durm(1)*100+durd(1),'%4.4d');
                tst2 = num2str(durm(30)*100+durd(30),'%4.4d');

                % xtlabels ={};
                % for tl = 1:30
                %     xtlabels{tl} =  [ num2str(durm(tl),'%2.2d') '/' num2str(durd(tl),'%2.2d') ];
                % end
                
                
                
                
                
                ffnm = ['/shared/dondo/home/ml423/world_mts_runs/Basins_events_chars/'];
                fid = fopen([ffnm 'Basin' num2str(gid) '_event' rainevent '.txt'],'r');
                if fid == -1
                    %fclose(fid);
                    continue
                end
                windowtimes = fscanf(fid,'%d\n',[4,1]);
                fclose(fid);
                if length(windowtimes)<4
                   continue 
                end
                %close all
                
                
                if 1==1
                    
                    %fh = APL_frontcut(bid,eid);
                    
                    tstr = [num2str(dury(1)) '-' num2str(durm(1),'%2.2d') '/' num2str(durd(1),'%2.2d') '-' num2str(durm(30),'%2.2d') '/' num2str(durd(30),'%2.2d')];
                    
                    
                    obstep = 90;
                    
                    tmp = load(['/shared/dondo/home/ml423/world_mts/obs_events/events/obs_' num2str(gid) '_event' rainevent '.mat']);
                    tmp = tmp.obsnow;
                    obsnow = tmp;
                    obsnan=find(isnan(obsnow));
                    if isempty(obsnan)
                    else
                        for na = 1:length(obsnan)
                            obsnow(obsnan(na)) = 0.5*(obsnow(obsnan(na)-1)+obsnow(obsnan(na)+1));
                        end
                    end
                    obs2 = interp1(1:31,obsnow,[1:0.33333:31]);
                    obsnow = obs2(1:90);
                    
                    % figure
                    % plot(obsnow)
                    
                    tmp = load(['/shared/dondo/home/ml423/world_mts/Precip/Event/ERA5_Basin' num2str(gid) '_' rainevent '.mat']);
                    s4dbkc = tmp.rain;

                    MCp = zeros(c,r);
                    pts = [];
                    dbkcr = 0;rmax = [];
                    for i = 1:720
                        % MCp = MCp + s4dbkc{i};
                        % pts(i) = mean(s4dbkc{i}(intt));
                        % dbkcr = dbkcr+pts(i);
                        rmax(i) = max(s4dbkc{i}(intt));
                    end
                    %MCp(ind) = nan;
                    %ori_mean = nanmean(MCp(:));
                    %ama = jet(40);ama(1,:) = [1 1 1];
                    % figure
                    % imagesc(MCp')
                    % colorbar
                    % colormap(ama)
                    % caxis([0.2*min(MCp(:)) 1.1*max(MCp(:))])
                    % 
                    
                    
                    % for Nature Science paper, consider the major events only
                    
                    launchp = find(rmax>0.5);% when rainfall >1mm/hr launch particles
                    ptmaxthre = 0.5*max(rmax);% used as a threshold to separate events
                    la2 = find(launchp<(windowtimes(4)*8));
                    la1 = find(launchp<(windowtimes(3)*8));
                    accu0 = 0;
                    if isempty(la1)
                        continue
                    end
                    for tmpla = launchp(la1(end)):-1:1
                        if ismember(tmpla,launchp)
                            accu0 = 0;
                            rainbegintime = tmpla;
                        else
                            accu0 = accu0+1;
                        end

                        if area>1000
                            noraininterval = 12;% for big basins, it is likely there is always rain in the basin, therefore using a small norain limit to distinguish different rain events.
                        else
                            noraininterval = 24;
                        end
                        if accu0>noraininterval
                            % consecutive x hours of basin max rain less than 0.5mm/hr,
                            % this should be counted as a different rain cell, should not
                            % be tracked!
                            % one exception is if another major rain is 48-24 ahead, should be
                            % considered as well
                            rsmallest = max([(tmpla-noraininterval),1]);

                            reconsider_period = rsmallest:tmpla;
                            %reconsider_period = (tmpla-noraininterval):tmpla;
                            secondary_rain = find(rmax(reconsider_period)>0.5);

                            if max(rmax(reconsider_period))>ptmaxthre
                                % adjust the rainbegintime here
                                rainbegintime = reconsider_period(secondary_rain(1));
                            end
                            smalltmpla = max([floor(tmpla/8),2]);
                            if obsnow(smalltmpla)>obsnow(smalltmpla-1)
                                % meaning obsflow is still rising.
                                if ~isempty(secondary_rain)
                                    rainbegintime = reconsider_period(secondary_rain(1));
                                end
                            end

                            break
                        end
                    end


                    
                    la0 = find(launchp == rainbegintime);

                    refinedlaunchp = launchp(la0:la2(end));
                    stat_left_bd_hydro = windowtimes(1);% systematic defined
                    stat_left_bd_rain = floor(rainbegintime/8);%target rainevent adjusted
                    stat_left_bd = max([stat_left_bd_hydro,stat_left_bd_rain]);

                    ICtiming = refinedlaunchp(1);%1-720, can be used to extract original IC from IRCICC studies

                    upcap = min([90 (windowtimes(4)+6)]);

                    if upcap<=stat_left_bd
                        continue
                    end

                    no3sum = nan(aiens,720);
                    %no2sum = [];
                    rts1 = eval_ind{bid,eid};
                    AI_mean_rain = zeros(c,r,length(rts1));
                    for replic = 1:aiens

                        ffnmx=['/taiga/ML/AI_20/Basin' num2str(gid) 'outputs/' rainevent '_AI_' num2str(replic,'%2.2d') '_' AItool '/' rainevent '/' flowtype{ftp}];
                        fid=fopen(ffnmx,'rb','ieee-le');
                        if fid == -1
                            disp('no AI done')
                            aibreak = 1;
                            break
                        else
                            aibreak = 0;
                        end
                        tmpdata=fread(fid,inf,'single');
                        fclose(fid);
                        if length(tmpdata)<(c*r*720)
                            disp([ii,ii])
                            unfinish_alps(ii,:) = [bid,eid];

                            aibreak = 1;
                            break
                        end
                        NODAc2o5=reshape(tmpdata,c,r,720);
                        no3 = [];
                        for k = 1:720
                            no3(k) = abs(NODAc2o5(cg,rg,k));
                        end

                        no3sum(replic,:) = no3;



                        % ffnm=['/taiga/ML/AI_20/Basin' num2str(gid) 'outputs/' rainevent '_AI_' num2str(replic,'%2.2d') '_' AItool '/' rainevent '/precipitation'];
                        % fid=fopen(ffnm,'rb','ieee-le');
                        % tmpdata=fread(fid,inf,'single');
                        % fclose(fid);
                        % NODAc2o5=reshape(tmpdata,c,r,720);
                        % flownoda5 = [];
                        % meanr250 = 0;
                        % tmpsum = zeros(c,r);
                        % for k = 1:720
                        %     tmp = NODAc2o5(:,:,k);
                        %     tmpsum = tmpsum+tmp;
                        %     flownoda5(k) = mean(tmp(intt));
                        %     meanr250 = meanr250+flownoda5(k);
                        % end
                        % no2 = flownoda5;%GREEN
                        % rainsd = std(tmpsum(intt));
                        % no2sum = [no2sum;no2];

                        % if replic == 1
                        %     figure
                        %     imagesc(tmpsum')
                        %     colorbar
                        %     colormap(jet(30))
                        % end
                        ffnmx=['/shared/dondo/home/ml423/world_mts/Basin' num2str(gid) '_input_ERA5_900m_tmp/' rainevent '/Precipitation_ERA5_lRFR_ori_' num2str(replic,'%3.3d') ];
                        fid=fopen(ffnmx,'rb','ieee-le');
                        tmpdata=fread(fid,inf,'single');
                        fclose(fid);
                        NODAc2o5=reshape(tmpdata,c,r,720);
                        
                        AI_mean_rain = AI_mean_rain+NODAc2o5(:,:,rts1);
                    end


                    if aibreak==1
                       continue 
                    end
                    rain_AI{bid,eid} = AI_mean_rain/aiens;
                    no3vol_rank = sum(no3sum,2);
                    [~,mloc] = sort(no3vol_rank);

                    if aiens>5
                        no3m = mean(no3sum(mloc(6:15),:),1);   
                    else
                        no3m = mean(no3sum,1);
                    end

                    

                    % filenam = str2num(rainevent)*100+fileind;
                    % ffnmx=['/shared/dondo/home/ml423/world_mts_runs/IRC1/Basin' num2str(gid) 'outputs/refxtmp' num2str(filenam) 'fullresults/' rainevent '_SW/' rainevent '/' flowtype{ftp}];
                    % fid=fopen(ffnmx,'rb','ieee-le');
                    % tmpdata=fread(fid,inf,'single');
                    % fclose(fid);
                    % NODAc2o5=reshape(tmpdata,c,r,720);
                    % flownoda5 = [];
                    % for k = 1:720
                    %     flownoda5(k) = abs(NODAc2o5(cg,rg,k));
                    % end
                    % no5 = flownoda5;% original era5 with original IC

                    ffnmx=['/taiga/ML/AI_20/Basin' num2str(gid) 'outputs/' rainevent '_AI_' num2str(0,'%2.2d') '/' rainevent '/' flowtype{ftp}];
                    fid=fopen(ffnmx,'rb','ieee-le');
                    tmpdata=fread(fid,inf,'single');
                    fclose(fid);
                    if length(tmpdata)<(c*r*720)
                       disp([ii,ii])
                    end
                    NODAc2o5=reshape(tmpdata,c,r,720);
                    flownoda5 = [];
                    for k = 1:720
                        flownoda5(k) = abs(NODAc2o5(cg,rg,k));
                    end
                    no6 = flownoda5;% original era5 with best IC from IRCicc


                    cal_window = (stat_left_bd-1)*8+1:(upcap-1)*8+1;
                    % mass balance check fo W00 and W34
                    %{
                        filenamp = str2num(rainevent)*100+fileind+2;% the very first P
                        cc = ls(['/shared/dondo/home/ml423/HMIOP/NoDAre/IRC_LB/Basin' num2str(bid,'%2.2d') 'outputs/refxtmp' num2str(filenamp) '/xtmpbb1*']);
                        cc2 = strsplit(cc);
                        tmpC  = struct2cell(load([cc2{end-1}]));
                        P00 = tmpC{1};
                        filenamp = str2num(rainevent)*100+fileind+24;% the very first P
                        cc = ls(['/shared/dondo/home/ml423/HMIOP/NoDAre/IRC_LB/Basin' num2str(bid,'%2.2d') 'outputs/refxtmp' num2str(filenamp) '/xtmpbb1*']);
                        cc2 = strsplit(cc);
                        tmpC  = struct2cell(load([cc2{1}]));
                        P34 = tmpC{1};
                        P00tot = 0;
                        P34tot = 0;
                        for k = 1:288
                            tmp = P00{k};
                            P00tot = P00tot+mean(tmp(intt));
                            tmp = P34{k};
                            P34tot = P34tot+mean(tmp(intt));
                        end
                        % mass balance quick calculations
                        runoff250 = sum(no3)*300/length(intt)/250/250*1000;



                        fixedstr = ['/shared/dondo/home/ml423/HMIOP/NoDAre/LB/APL_Basin' num2str(bid,'%2.2d') '_fixed_250m/'];
                        datastr = ['/shared/dondo/home/ml423/HMIOP/NoDAre/IRC_LB/Basin' num2str(bid,'%2.2d') 'outputs/refxtmp' num2str(filenam) 'fullresults/Basin' num2str(bid,'%2.2d') '_output_250m5min_' rainevent '_SW/' rainevent '00/'];

                        lastnote = 'c';
                        if p==0
                            laststr = ['/shared/dondo/home/ml423/HMIOP/NoDAre/LB/Basin' num2str(bid,'%2.2d') '_output_250m5min_' rainevent '_' lastnote '/'];
                            [smdelta250,gddelta250] = waterbalance(0,fixedstr,c,r,288,288,laststr,datastr,'',intt,250);
                        elseif p==9
                            if ICC==1% if ICC is part of the algo
                                laststr = ['/shared/dondo/home/ml423/HMIOP/NoDAre/IRC_LB/ICall/' num2str(ICn) '/Basin' num2str(bid,'%2.2d') '_output_250m5min_' rainevent lastnote '/'];
                            else
                                laststr = ['/shared/dondo/home/ml423/HMIOP/NoDAre/LB/Basin' num2str(bid,'%2.2d') '_output_250m5min_' rainevent '_' lastnote '/'];
                            end
                            [smdelta250,gddelta250] = waterbalance(0,fixedstr,c,r,288,288,laststr,datastr,'',intt,250);
                        end

                        if p==0
                            rainDBKC250 = dbkcr;% in mm
                            riverDBKC250 = runoff250;% in mm
                            deltasDBKC250 = smdelta250/250/250/length(intt)*1000; %in mm
                            gdDBKC250 = gddelta250/250/250/length(intt)*1000; %in mm
                            budgetDBKC250 = rainDBKC250-riverDBKC250-(deltasDBKC250+gdDBKC250);
                            wbterms_DBKC = [rainDBKC250,riverDBKC250,deltasDBKC250,gdDBKC250,budgetDBKC250];
                            budgetalldbkc = [budgetalldbkc;[str2num(rainevent),wbterms_DBKC]];
                        elseif p==9
                            rainIRC250 = P34tot;% in mm
                            riverIRC250 = runoff250;% in mm
                            deltasIRC250 = smdelta250/250/250/length(intt)*1000; %in mm
                            gdIRC250 = gddelta250/250/250/length(intt)*1000; %in mm
                            budgetIRC250 = rainIRC250-riverIRC250-(deltasIRC250+gdIRC250);
                            wbterms_IRC = [rainIRC250,riverIRC250,deltasIRC250,gdIRC250,budgetIRC250];
                            budgetallirc = [budgetallirc;[str2num(rainevent),wbterms_IRC]];
                        end


                    %}


                    simu = no6(1:8:end);


                    simustat = simu(stat_left_bd:upcap);
                    obsnowstat = obsnow(stat_left_bd:upcap);
                    tmp = 0;
                    tmp1 = 0;
                    tmp2 = mean(obsnowstat);
                    for i = 1:length(obsnowstat)
                        tmp = tmp + (simustat(i)-obsnowstat(i))^2;
                        tmp1 = tmp1 + (tmp2-obsnowstat(i))^2;
                    end
                    NSEori = 1-tmp./tmp1;

                    rr = corrcoef(simustat',obsnowstat);
                    rstar = rr(1,2);
                    obsstd = std(obsnowstat);
                    simustd = std(simustat);
                    obsmean = mean(obsnowstat);
                    simumean = mean(simustat);
                    KGEori = 1-sqrt((rstar-1)^2+(simustd/obsstd-1)^2+(simumean/obsmean-1)^2);

                    simu = no3m(1:8:end);
                    simustat = simu(stat_left_bd:upcap);
                    obsnowstat = obsnow(stat_left_bd:upcap);
                    tmp = 0;
                    tmp1 = 0;
                    tmp2 = mean(obsnowstat);
                    for i = 1:length(obsnowstat)
                        tmp = tmp + (simustat(i)-obsnowstat(i))^2;
                        tmp1 = tmp1 + (tmp2-obsnowstat(i))^2;
                    end
                    NSE = 1-tmp./tmp1;

                    rr = corrcoef(simustat',obsnowstat);
                    rstar = rr(1,2);
                    obsstd = std(obsnowstat);
                    simustd = std(simustat);
                    obsmean = mean(obsnowstat);
                    simumean = mean(simustat);
                    KGE = 1-sqrt((rstar-1)^2+(simustd/obsstd-1)^2+(simumean/obsmean-1)^2);

                    [vsimu,tsimu]=max(simustat);
                    [vobs,tobs]=max(obsnowstat);
                    EPV = vsimu-vobs;
                    EPT = (tsimu-tobs)*1;
                    simvol = 0;
                    obsvol = 0;
                    for vlo=2:length(obsnowstat)-1
                        obsvol = obsvol+obsnowstat(vlo)*1;
                        simvol = simvol+simustat(vlo)*1;
                    end
                    obsvol = obsvol+0.5*(obsnowstat(1)+obsnowstat(end))*1;
                    simvol = simvol+0.5*(simustat(1)+simustat(end))*1;

                    EV = (simvol-obsvol)./obsvol;
                    alps_AIstats(ii,:) = [bid,eid,KGEori,KGE];
                 
                    
                    %{
                
                    f = figure('visible','on');
                    f.Position = [5 5 700 350];
%                     
                    pos1 = [0.13 0.8 0.67 0.15];
                    %subplot('Position',pos1)
                    pos2 = [0.15 0.32 0.78 0.42];
                    
                    subplot('Position',pos2)
                    % hjh = plot(2:3:obstep*3,obsdatav1(:,6),'k',1:288,no3,'g',1:288,no2,'r',1:288,no5,'m');
                    % set(hjh,'LineWidth',2.5)
                    % legend('Obs.','IMERG_D','STIV_D','STIV_{DBK}','Location','Northeast')
                    %yyaxis left
                    xj = 2:8:obstep*8;
                    xj = xj';
                    bd1 = 1.05*obsnow;
                    bd2 = 0.95*obsnow;
                    sh = scatter(2:8:obstep*8,obsnow,45,'g','filled');
                    sh.MarkerEdgeColor = [0,1,0];
                    for ii = 1:size(no3sum,1)
                        hold on
                        hjh = plot(1:8:720,no3sum(ii,1:8:720),'-.','LineWidth',1.5,'Color',[0, 0, 1, 0.3]);
                    end
                   
%                     hold on
%                     hjh = plot(1:8:720,no5(1:8:720),'m:','LineWidth',2);
                    hold on
                    hjh = plot(1:8:720,no6(1:8:720),'r:','LineWidth',2);
                    hold on
                    hjh = plot(1:8:720,no3m(1:8:720),'b-.');
                    hold on
                    hjb = plot([1+8*stat_left_bd 1+8*stat_left_bd],[0 1000000],'k-.');
                    hold on
                    hjb2 = plot([1+8*(windowtimes(4)+6) 1+8*(windowtimes(4)+6)],[0 1000000],'k-.');
                    
                    % ymax is based on max of y but rounded to
                    % nearest interval defined by max of y
                   
                    ymm = ceil(max(max(1.5*obsnow),max(1.5*no3m(:)))/(max(obsnow)/4))*(max(obsnow)/4);
                    
                    tf = text(16,0.85*ymm,['' num2str(KGEori,'%2.2f')],'Color','r');
                    tf.FontSize = 14;
                    tf.FontWeight = 'bold';
                    
                    tf = text(16,0.72*ymm,['' num2str(KGE,'%2.2f')],'Color','b');
                    tf.FontSize = 14;
                    tf.FontWeight = 'bold';
                    set(hjh,'LineWidth',2.5)
                    %             hold on
                    %             patch([xj; flipud(xj)]', [bd1; flipud(bd2)]',  'k', 'FaceAlpha',0.2, 'EdgeColor','none')
                    %             hold off
                    ylim([0 ymm])
                    ylabel('Streamflow (m^3/s)')
                    
                    box on
                    
                    %title([num2str(gid) ' IRC W' num2str(itranum) num2str(wnum) ' | ' num2str(area,'%0.0f') ' km^2'])
                    
                    
                    xlabel([tstr ' '])
                    %set(gca,'XTick',1:3:288,'XTickLabel',[])
                    xlim([0 719])
                    set(gca,'XTick',0:144:720,'XTickLabel',xtlabels(1:6:30),'FontSize',12,'LineWidth',2.5,'FontWeight','bold')
                    
                    
                    hax = gca;
                    hax.XAxis.MinorTickValues = linspace(0,720,31);
                    hax.XMinorTick = 'on';
                    
                    
                    set(gca,'FontName','Times New Roman','FontSize',14,'LineWidth',2.5,'FontWeight','bold');
                    Tight = get(gca, 'TightInset');  %Gives you the bording spacing between plot box and any axis labels
                    %[Left Bottom Right Top] spacing
                    NewPos = [Tight(1) Tight(2) 1-Tight(1)-Tight(3)-0.01 1-Tight(2)-Tight(4)-0.01]; %New plot position [X Y W H]
                    set(gca, 'Position', NewPos);
                    %saveas(gcf,['/shared/dondo/home/ml423/world_mts_runs/figs/IRC_images/Basin' num2str(gid) '_' rainevent '_' num2str(ICn) '_W' num2str(itranum) num2str(wnum) '.jpg'])
                    %exportgraphics(gcf,['/taiga/ML/AI_20/pics/Basin' num2str(gid) '_ind_' num2str(bid) 'event' num2str(eid,'%2.2d') '_AI_20ens_' AItool '_20.jpg'],'Resolution',300)

                    %}
                end
                
                %}  
            end
            % done plot optimum Wxy
            
            
        end
    end
end
% 

ubid = unique(alps_AIstats(:,1));ubid(find(ubid==0)) = [];
ubidinfo = [];ubidfull = {};
for i = 1:length(ubid)
    tmp = find(alps_AIstats(:,1)==ubid(i));
    AIori = nanmedian(alps_AIstats(tmp,3));
    AIopt = nanmedian(alps_AIstats(tmp,4));
    goodratio = length(find(alps_AIstats(tmp,4)>alps_AIstats(tmp,3)))/length(tmp);
    ubidinfo(i,:) = [ubid(i),AIori,AIopt,goodratio];
    ubidfull{i} = [ones(length(tmp),1)*ubid(i),alps_AIstats(tmp,2),WORLD_events(ubid(i),1+alps_AIstats(tmp,2))',alps_AIstats(tmp,3),alps_AIstats(tmp,4)];
    % figure
    % cdfplot(ubidfull{i}(:,2))
    % hold on
    % cdfplot(ubidfull{i}(:,3))
    % title(num2str(ubidfull{i}(1)))
end

%% plot AI outputs for Andes 

clear tmpp
clear tmp
clear predatatmp
clear data3D
addpath('/shared/dondo/home/ml423/HMIOP/NoDAre/New era/')
addpath('/shared/dondo/home/ml423/DA/Critical matrices/')
addpath('/shared/dondo/home/ml423/DA/slanCM/')

tmp = load(['/shared/dondo/home/ml423/DA/Critical matrices/world_selected_gauge_events_allinfo_v3.mat'],'gauge_events_first_last');
Basins_sum = tmp.gauge_events_first_last;clear tmp;
tmp = load(['/shared/dondo/home/ml423/DA/Critical matrices/WORLD_events_v3.mat'],'WORLD_events');
WORLD_events = tmp.WORLD_events;clear tmp;

flowtype = {'streamflowarea','interflowarea','overlandflowarea','baseflowarea'};
ftp = 1;
AItool = 'A15';

%hima_AIstats = [];

% figure
% cdfplot(hima_AIstats(831:880,3))
% hold on
% cdfplot(hima_AIstats(831:880,4))
% figure
% scatter(hima_AIstats(:,3),hima_AIstats(:,4),18,hima_AIstats(:,1),'filled')
% caxis([182 250])
% colormap(distinguishable_colors(30))
% xlim([-2 1])
% ylim([-2 1])
andessmall = find(gauge_events_first_last(:,5)<1200);
adsall = [];
tooflat = [215 242 245 296 330 368 390];%critical
for ki = 1:length(andessmall)
    bid = find(WORLD_events(:,1)==gauge_events_first_last(andessmall(ki),1));
    if bid>=134&bid<=411
        if ~ismember(bid,tooflat)
            if gauge_events_first_last(andessmall(ki),3)>19900000
                adsall = [adsall;bid];
            end
        end
    end
end
andes_alleve = [];
for i = 1:length(adsall)
    bid = adsall(i);
    ixx = find(~isnan(WORLD_events(bid,:)));
    etot = length(ixx)-1;
    if etot == 0
        continue
    else
        if etot>149
            etot = 149
        end
        for eid = 1:etot
            close all
            % the if comments below are data quality control of the Alps
            % and hima
            %if isnan(event_quality(bid,eid+1))||isnan(kge_optm(bid,eid))
            if eqcontrol(bid,eid)==1%isempty(dssumtmp{bid,eid})||cvtmp(bid,eid)>1000
                andes_alleve = [andes_alleve;[bid,eid]];
            end
        end
    end
end


aiens = 5;
andes_AIstats = [];
clear unfinish_andes
for ii = 1:length(andes_alleve)%[329 332 333 334 335 336 337 338 184 185 186 187 290 296 298 299 301 302 215 216 217 224 229 241 242 250 266 299 320 321 322 314]%[215 184]%[3 50]
    
    bid = andes_alleve(ii,1)
    gid = WORLD_events(bid,1);
    
    gloc = find(Basins_sum(:,1)==gid);
    
    r = Basins_sum(gloc,6);
    c = Basins_sum(gloc,7);
    rg = Basins_sum(gloc,8);
    cg = Basins_sum(gloc,9);
    area = Basins_sum(gloc,5);
    
    % if bid>=412
    %     cvtmp = cvfloodhima;
    %     dssumtmp = dssumhima;
    % elseif bid>=134
    %     cvtmp = cvfloodandes;
    %     dssumtmp = dssumandes;
    % end
    
    tmp = load(['/shared/dondo/home/ml423/DA/Critical matrices/World_Basins/inttind/Basin' num2str(gid) '_ind.mat'],'ind');
    ind = tmp.ind;
    tmp = load(['/shared/dondo/home/ml423/DA/Critical matrices/World_Basins/inttind/Basin' num2str(gid) '_intt.mat'],'intt');
    intt = tmp.intt;
    % ffnm=['/shared/dondo/home/ml423/world_mts/Basin' num2str(gid) '_fixed_900m/dem.bin'];
    % fid=fopen(ffnm,'rb','ieee-le');
    % tmpdata=fread(fid,inf,'single');
    % fclose(fid);
    % dem = reshape(tmpdata,c,r);
    % dem(ind) = nan;
    % demf = dem';
    % colo = slanCM('terrain',80);
    % colo(1,:) = [0.8 0.8 0.8];
    % figure
    % imagesc(demf)
    % colormap(colo)
    
    % ffnm=['/shared/dondo/home/ml423/world_mts/Basin' num2str(gid) '_fixed_900m/countfile.out'];
    % fid=fopen(ffnm);
    % facc=fscanf(fid,'%d',[c,r]);
    % fclose(fid);
    % facc(ind) = nan;
    % [streamx,streamy] = find(facc>5);
    
    
    
    ICn = 4140;
    ICC = 0;% IF ICC is part of the algo
    
    fileind = -100000*ICn;
    KGEsu = [];
    evesu = [];
    budgetalldbkc = [];
    budgetallirc = [];
    
    ixx = find(~isnan(WORLD_events(bid,:)));
    etot = length(ixx)-1;
    if etot == 0
        continue
    else
        if etot>149
           etot = 149; 
        end
        for eidd = 1
            eid = andes_alleve(ii,2);
            %close all
            % the if comments below are data quality control of the Alps
            % and hima
            %if isnan(event_quality(bid,eid+1))||isnan(kge_optm(bid,eid))
            if eqcontrol(bid,eid)~=1%isempty(dssumtmp{bid,eid})||cvtmp(bid,eid)>1000    
                disp('Event quality is bad, peak too early or too late')
                continue
            else
                rainevent = num2str(WORLD_events(bid,eid+1));
                yth = floor(WORLD_events(bid,eid+1)/10000);
                dth = mod(WORLD_events(bid,eid+1),100);
                mth = floor((WORLD_events(bid,eid+1)-yth*10000)/100);
                
                thresdate = datetime(yth,mth,dth);
                
                durati = (thresdate-days(15)):days(1):(thresdate+days(15));
                [dury,durm,durd] = ymd(durati);
                tst1 = num2str(durm(1)*100+durd(1),'%4.4d');
                tst2 = num2str(durm(30)*100+durd(30),'%4.4d');
                
                xtlabels ={};
                for tl = 1:30
                    xtlabels{tl} =  [ num2str(durm(tl),'%2.2d') '/' num2str(durd(tl),'%2.2d') ];
                end
                
                
                
                
                
                ffnm = ['/shared/dondo/home/ml423/world_mts_runs/Basins_events_chars/'];
                fid = fopen([ffnm 'Basin' num2str(gid) '_event' rainevent '.txt'],'r');
                if fid == -1
                    %fclose(fid);
                    continue
                end
                windowtimes = fscanf(fid,'%d\n',[4,1]);
                fclose(fid);
                if length(windowtimes)<4
                   continue 
                end
                %close all
                
                
                if 1==1
                    
                    %fh = APL_frontcut(bid,eid);
                    
                    tstr = [num2str(dury(1)) '-' num2str(durm(1),'%2.2d') '/' num2str(durd(1),'%2.2d') '-' num2str(durm(30),'%2.2d') '/' num2str(durd(30),'%2.2d')];
                    
                    
                    obstep = 90;
                    
                    tmp = load(['/shared/dondo/home/ml423/world_mts/obs_events/events/obs_' num2str(gid) '_event' rainevent '.mat']);
                    tmp = tmp.obsnow;
                    obsnow = tmp;
                    obsnan=find(isnan(obsnow));
                    if isempty(obsnan)
                    else
                        for na = 1:length(obsnan)
                            obsnow(obsnan(na)) = 0.5*(obsnow(obsnan(na)-1)+obsnow(obsnan(na)+1));
                        end
                    end
                    obs2 = interp1(1:31,obsnow,[1:0.33333:31]);
                    obsnow = obs2(1:90);
                    
                    %figure
                    %plot(obsnow)
                    
                    tmp = load(['/shared/dondo/home/ml423/world_mts/Precip/Event/ERA5_Basin' num2str(gid) '_' rainevent '.mat']);
                    s4dbkc = tmp.rain;
                    
                    MCp = zeros(c,r);
                    pts = [];
                    dbkcr = 0;rmax = [];
                    for i = 1:720
                        % MCp = MCp + s4dbkc{i};
                        % pts(i) = mean(s4dbkc{i}(intt));
                        % dbkcr = dbkcr+pts(i);
                        rmax(i) = max(s4dbkc{i}(intt));
                    end
                    % MCp(ind) = nan;
                    % ori_mean = nanmean(MCp(:));
                    % ama = jet(40);ama(1,:) = [1 1 1];
                    % figure
                    % imagesc(MCp')
                    % colorbar
                    % colormap(ama)
                    % caxis([0.2*min(MCp(:)) 1.1*max(MCp(:))])
                    % 
                    
                    
                    % for Nature Science paper, consider the major events only
                    
                    launchp = find(rmax>0.5);% when rainfall >1mm/hr launch particles
                    ptmaxthre = 0.5*max(rmax);% used as a threshold to separate events
                    la2 = find(launchp<(windowtimes(4)*8));
                    la1 = find(launchp<(windowtimes(3)*8));
                    accu0 = 0;
                    if isempty(la1)
                        continue
                    end
                    for tmpla = launchp(la1(end)):-1:1
                        if ismember(tmpla,launchp)
                            accu0 = 0;
                            rainbegintime = tmpla;
                        else
                            accu0 = accu0+1;
                        end
                        
                        if area>1000
                            noraininterval = 12;% for big basins, it is likely there is always rain in the basin, therefore using a small norain limit to distinguish different rain events.
                        else
                            noraininterval = 24;
                        end
                        if accu0>noraininterval
                            % consecutive x hours of basin max rain less than 0.5mm/hr,
                            % this should be counted as a different rain cell, should not
                            % be tracked!
                            % one exception is if another major rain is 48-24 ahead, should be
                            % considered as well
                            rsmallest = max([(tmpla-noraininterval),1]);
                            
                            reconsider_period = rsmallest:tmpla;
                            %reconsider_period = (tmpla-noraininterval):tmpla;
                            secondary_rain = find(rmax(reconsider_period)>0.5);
                            
                            if max(rmax(reconsider_period))>ptmaxthre
                                % adjust the rainbegintime here
                                rainbegintime = reconsider_period(secondary_rain(1));
                            end
                            smalltmpla = max([floor(tmpla/8),2]);
                            if obsnow(smalltmpla)>obsnow(smalltmpla-1)
                                % meaning obsflow is still rising.
                                if ~isempty(secondary_rain)
                                    rainbegintime = reconsider_period(secondary_rain(1));
                                end
                            end
                            
                            break
                        end
                    end
                    
                    la0 = find(launchp == rainbegintime);
                    
                    refinedlaunchp = launchp(la0:la2(end));
                    stat_left_bd_hydro = windowtimes(1);% systematic defined
                    stat_left_bd_rain = floor(rainbegintime/8);%target rainevent adjusted
                    stat_left_bd = max([stat_left_bd_hydro,stat_left_bd_rain]);
                    
                    ICtiming = refinedlaunchp(1);%1-720, can be used to extract original IC from IRCICC studies
                    
                    upcap = min([90 (windowtimes(4)+6)]);
                    
                    if upcap<=stat_left_bd
                        continue
                    end
                    
                    no3sum = nan(aiens,720);
                    %no2sum = [];
                    rts1 = eval_ind{bid,eid};
                    AI_mean_rain = zeros(c,r,length(rts1));
                    for replic = 1:aiens
                        
                        ffnmx=['/taiga/ML/AI_20/Basin' num2str(gid) 'outputs/' rainevent '_AI_' num2str(replic,'%2.2d') '_' AItool '/' rainevent '/' flowtype{ftp}];
                        fid=fopen(ffnmx,'rb','ieee-le');
                        if fid == -1
                            disp('no AI done')
                            aibreak = 1;
                            break
                        else
                            aibreak = 0;
                        end
                        tmpdata=fread(fid,inf,'single');
                        fclose(fid);
                        if length(tmpdata)<(c*r*720)
                            disp([ii,ii])
                            unfinish_andes(ii,:) = [bid,eid];
                            aibreak = 1;
                            break
                        end
                        NODAc2o5=reshape(tmpdata,c,r,720);
                        no3 = [];
                        for k = 1:720
                            no3(k) = abs(NODAc2o5(cg,rg,k));
                        end
                       
                        no3sum(replic,:) = no3;
                        % ffnm=['/taiga/ML/AI_20/Basin' num2str(gid) 'outputs/' rainevent '_AI_' num2str(replic,'%2.2d') '_' AItool '/' rainevent '/precipitation'];
                        % fid=fopen(ffnm,'rb','ieee-le');
                        % tmpdata=fread(fid,inf,'single');
                        % fclose(fid);
                        % NODAc2o5=reshape(tmpdata,c,r,720);
                        % flownoda5 = [];
                        % meanr250 = 0;
                        % tmpsum = zeros(c,r);
                        % for k = 1:720
                        %     tmp = NODAc2o5(:,:,k);
                        %     tmpsum = tmpsum+tmp;
                        %     flownoda5(k) = mean(tmp(intt));
                        %     meanr250 = meanr250+flownoda5(k);
                        % end
                        % no2 = flownoda5;%GREEN
                        % rainsd = std(tmpsum(intt));
                        % no2sum = [no2sum;no2];
                        % 
                        % if replic == 1
                        %     figure
                        %     imagesc(tmpsum')
                        %     colorbar
                        %     colormap(jet(30))
                        % end
                        ffnmx=['/shared/dondo/home/ml423/world_mts/Basin' num2str(gid) '_input_ERA5_900m_tmp/' rainevent '/Precipitation_ERA5_lRFR_ori_' num2str(replic,'%3.3d') ];
                        fid=fopen(ffnmx,'rb','ieee-le');
                        tmpdata=fread(fid,inf,'single');
                        fclose(fid);
                        NODAc2o5=reshape(tmpdata,c,r,720);
                        
                        AI_mean_rain = AI_mean_rain+NODAc2o5(:,:,rts1);

                    end
                    
                    if aibreak==1
                       continue 
                    end
                    rain_AI{bid,eid} = AI_mean_rain/aiens;
                    no3vol_rank = sum(no3sum,2);
                    [~,mloc] = sort(no3vol_rank);
                    if aiens>5
                       no3m = mean(no3sum(mloc(6:15),:),1);
                    else
                       no3m = mean(no3sum,1); 
                    end

                    
%                     filenam = str2num(rainevent)*100+fileind;
%                     ffnmx=['/shared/dondo/home/ml423/world_mts_runs/IRC1/Basin' num2str(gid) 'outputs/refxtmp' num2str(filenam) 'fullresults/' rainevent '_SW/' rainevent '/' flowtype{ftp}];
%                     fid=fopen(ffnmx,'rb','ieee-le');
%                     tmpdata=fread(fid,inf,'single');
%                     fclose(fid);
%                     NODAc2o5=reshape(tmpdata,c,r,720);
%                     flownoda5 = [];
%                     for k = 1:720
%                         flownoda5(k) = abs(NODAc2o5(cg,rg,k));
%                     end
%                     no5 = flownoda5;% original era5 with original IC
%                     
                    ffnmx=['/taiga/ML/AI_20/Basin' num2str(gid) 'outputs/' rainevent '_AI_' num2str(0,'%2.2d') '/' rainevent '/' flowtype{ftp}];
                    fid=fopen(ffnmx,'rb','ieee-le');
                    % if fid == -1
                    %     continue
                    % end
                    tmpdata=fread(fid,inf,'single');
                    fclose(fid);
                    NODAc2o5=reshape(tmpdata,c,r,720);
                    flownoda5 = [];
                    for k = 1:720
                        flownoda5(k) = abs(NODAc2o5(cg,rg,k));
                    end
                    no6 = flownoda5;% original era5 with best IC from IRCicc
                    
                    
                    cal_window = (stat_left_bd-1)*8+1:(upcap-1)*8+1;
                    % mass balance check fo W00 and W34
                    %{
                        filenamp = str2num(rainevent)*100+fileind+2;% the very first P
                        cc = ls(['/shared/dondo/home/ml423/HMIOP/NoDAre/IRC_LB/Basin' num2str(bid,'%2.2d') 'outputs/refxtmp' num2str(filenamp) '/xtmpbb1*']);
                        cc2 = strsplit(cc);
                        tmpC  = struct2cell(load([cc2{end-1}]));
                        P00 = tmpC{1};
                        filenamp = str2num(rainevent)*100+fileind+24;% the very first P
                        cc = ls(['/shared/dondo/home/ml423/HMIOP/NoDAre/IRC_LB/Basin' num2str(bid,'%2.2d') 'outputs/refxtmp' num2str(filenamp) '/xtmpbb1*']);
                        cc2 = strsplit(cc);
                        tmpC  = struct2cell(load([cc2{1}]));
                        P34 = tmpC{1};
                        P00tot = 0;
                        P34tot = 0;
                        for k = 1:288
                            tmp = P00{k};
                            P00tot = P00tot+mean(tmp(intt));
                            tmp = P34{k};
                            P34tot = P34tot+mean(tmp(intt));
                        end
                        % mass balance quick calculations
                        runoff250 = sum(no3)*300/length(intt)/250/250*1000;
                        
                        
                        
                        fixedstr = ['/shared/dondo/home/ml423/HMIOP/NoDAre/LB/APL_Basin' num2str(bid,'%2.2d') '_fixed_250m/'];
                        datastr = ['/shared/dondo/home/ml423/HMIOP/NoDAre/IRC_LB/Basin' num2str(bid,'%2.2d') 'outputs/refxtmp' num2str(filenam) 'fullresults/Basin' num2str(bid,'%2.2d') '_output_250m5min_' rainevent '_SW/' rainevent '00/'];
                        
                        lastnote = 'c';
                        if p==0
                            laststr = ['/shared/dondo/home/ml423/HMIOP/NoDAre/LB/Basin' num2str(bid,'%2.2d') '_output_250m5min_' rainevent '_' lastnote '/'];
                            [smdelta250,gddelta250] = waterbalance(0,fixedstr,c,r,288,288,laststr,datastr,'',intt,250);
                        elseif p==9
                            if ICC==1% if ICC is part of the algo
                                laststr = ['/shared/dondo/home/ml423/HMIOP/NoDAre/IRC_LB/ICall/' num2str(ICn) '/Basin' num2str(bid,'%2.2d') '_output_250m5min_' rainevent lastnote '/'];
                            else
                                laststr = ['/shared/dondo/home/ml423/HMIOP/NoDAre/LB/Basin' num2str(bid,'%2.2d') '_output_250m5min_' rainevent '_' lastnote '/'];
                            end
                            [smdelta250,gddelta250] = waterbalance(0,fixedstr,c,r,288,288,laststr,datastr,'',intt,250);
                        end
                        
                        if p==0
                            rainDBKC250 = dbkcr;% in mm
                            riverDBKC250 = runoff250;% in mm
                            deltasDBKC250 = smdelta250/250/250/length(intt)*1000; %in mm
                            gdDBKC250 = gddelta250/250/250/length(intt)*1000; %in mm
                            budgetDBKC250 = rainDBKC250-riverDBKC250-(deltasDBKC250+gdDBKC250);
                            wbterms_DBKC = [rainDBKC250,riverDBKC250,deltasDBKC250,gdDBKC250,budgetDBKC250];
                            budgetalldbkc = [budgetalldbkc;[str2num(rainevent),wbterms_DBKC]];
                        elseif p==9
                            rainIRC250 = P34tot;% in mm
                            riverIRC250 = runoff250;% in mm
                            deltasIRC250 = smdelta250/250/250/length(intt)*1000; %in mm
                            gdIRC250 = gddelta250/250/250/length(intt)*1000; %in mm
                            budgetIRC250 = rainIRC250-riverIRC250-(deltasIRC250+gdIRC250);
                            wbterms_IRC = [rainIRC250,riverIRC250,deltasIRC250,gdIRC250,budgetIRC250];
                            budgetallirc = [budgetallirc;[str2num(rainevent),wbterms_IRC]];
                        end
                        
                        
                    %}
                    
                    
                    simu = no6(1:8:end);
                    
                    
                    simustat = simu(stat_left_bd:upcap);
                    obsnowstat = obsnow(stat_left_bd:upcap);
                    tmp = 0;
                    tmp1 = 0;
                    tmp2 = mean(obsnowstat);
                    for i = 1:length(obsnowstat)
                        tmp = tmp + (simustat(i)-obsnowstat(i))^2;
                        tmp1 = tmp1 + (tmp2-obsnowstat(i))^2;
                    end
                    NSEori = 1-tmp./tmp1;
                    
                    rr = corrcoef(simustat',obsnowstat);
                    rstar = rr(1,2);
                    obsstd = std(obsnowstat);
                    simustd = std(simustat);
                    obsmean = mean(obsnowstat);
                    simumean = mean(simustat);
                    KGEori = 1-sqrt((rstar-1)^2+(simustd/obsstd-1)^2+(simumean/obsmean-1)^2);
                    
                    simu = no3m(1:8:end);
                    simustat = simu(stat_left_bd:upcap);
                    obsnowstat = obsnow(stat_left_bd:upcap);
                    tmp = 0;
                    tmp1 = 0;
                    tmp2 = mean(obsnowstat);
                    for i = 1:length(obsnowstat)
                        tmp = tmp + (simustat(i)-obsnowstat(i))^2;
                        tmp1 = tmp1 + (tmp2-obsnowstat(i))^2;
                    end
                    NSE = 1-tmp./tmp1;
                    
                    rr = corrcoef(simustat',obsnowstat);
                    rstar = rr(1,2);
                    obsstd = std(obsnowstat);
                    simustd = std(simustat);
                    obsmean = mean(obsnowstat);
                    simumean = mean(simustat);
                    KGE = 1-sqrt((rstar-1)^2+(simustd/obsstd-1)^2+(simumean/obsmean-1)^2);
              
                    [vsimu,tsimu]=max(simustat);
                    [vobs,tobs]=max(obsnowstat);
                    EPV = vsimu-vobs;
                    EPT = (tsimu-tobs)*1;
                    simvol = 0;
                    obsvol = 0;
                    for vlo=2:length(obsnowstat)-1
                        obsvol = obsvol+obsnowstat(vlo)*1;
                        simvol = simvol+simustat(vlo)*1;
                    end
                    obsvol = obsvol+0.5*(obsnowstat(1)+obsnowstat(end))*1;
                    simvol = simvol+0.5*(simustat(1)+simustat(end))*1;
                    
                    EV = (simvol-obsvol)./obsvol;
                    andes_AIstats(ii,:) = [bid,eid,KGEori,KGE];
                    %{
                    % KGE is period defined, while ther other stats
                    % use the whole series
                    
                    %                             disp(['NSE is ' num2str(NSE)])
                    %                             disp(['KGE is ' num2str(KGE)])
                    %                             disp(['EPV is ' num2str(EPV)])
                    %                             disp(['EPT is ' num2str(EPT)])
                    
                    
%                     
%                                                 figure
%                                                 t = tiledlayout(1,1,'Padding','none');
%                                                 t.Units = 'inches';
%                                                 t.OuterPosition = [0.25 0.25 4.8 3.2];
%                                                 nexttile;
                    
                    f = figure('visible','on');
                    f.Position = [5 5 700 350];
%                     
                    pos1 = [0.13 0.8 0.67 0.15];
                    %subplot('Position',pos1)
                    pos2 = [0.15 0.32 0.78 0.42];
                    
                    subplot('Position',pos2)
                    % hjh = plot(2:3:obstep*3,obsdatav1(:,6),'k',1:288,no3,'g',1:288,no2,'r',1:288,no5,'m');
                    % set(hjh,'LineWidth',2.5)
                    % legend('Obs.','IMERG_D','STIV_D','STIV_{DBK}','Location','Northeast')
                    %yyaxis left
                    xj = 2:8:obstep*8;
                    xj = xj';
                    bd1 = 1.05*obsnow;
                    bd2 = 0.95*obsnow;
                    sh = scatter(2:8:obstep*8,obsnow,45,'g','filled');
                    sh.MarkerEdgeColor = [0,1,0];
                    for ii = 1:size(no3sum,1)
                        hold on
                        hjh = plot(1:8:720,no3sum(ii,1:8:720),'-.','LineWidth',1.5,'Color',[0, 0, 1, 0.3]);
                    end
                   
%                     hold on
%                     hjh = plot(1:8:720,no5(1:8:720),'m:','LineWidth',2);
                    hold on
                    hjh = plot(1:8:720,no6(1:8:720),'r:','LineWidth',2);
                    hold on
                    hjh = plot(1:8:720,no3m(1:8:720),'b-.');
                    hold on
                    hjb = plot([1+8*stat_left_bd 1+8*stat_left_bd],[0 1000000],'k-.');
                    hold on
                    hjb2 = plot([1+8*(windowtimes(4)+6) 1+8*(windowtimes(4)+6)],[0 1000000],'k-.');
                    
                    % ymax is based on max of y but rounded to
                    % nearest interval defined by max of y
                   
                        ymm = ceil(max(max(1.5*obsnow),max(1.5*no3m(:)))/(max(obsnow)/4))*(max(obsnow)/4);
                    
                    tf = text(16,0.85*ymm,['' num2str(KGEori,'%2.2f')],'Color','r');
                    tf.FontSize = 14;
                    tf.FontWeight = 'bold';
                    
                    tf = text(16,0.72*ymm,['' num2str(KGE,'%2.2f')],'Color','b');
                    tf.FontSize = 14;
                    tf.FontWeight = 'bold';
                    set(hjh,'LineWidth',2.5)
                    %             hold on
                    %             patch([xj; flipud(xj)]', [bd1; flipud(bd2)]',  'k', 'FaceAlpha',0.2, 'EdgeColor','none')
                    %             hold off
                    ylim([0 ymm])
                    ylabel('Streamflow (m^3/s)')
                    
                    box on
                    
                    %title([num2str(gid) ' IRC W' num2str(itranum) num2str(wnum) ' | ' num2str(area,'%0.0f') ' km^2'])
                    
                    
                    xlabel([tstr ' '])
                    %set(gca,'XTick',1:3:288,'XTickLabel',[])
                    xlim([0 719])
                    set(gca,'XTick',0:144:720,'XTickLabel',xtlabels(1:6:30),'FontSize',12,'LineWidth',2.5,'FontWeight','bold')
                    
                    
                    hax = gca;
                    hax.XAxis.MinorTickValues = linspace(0,720,31);
                    hax.XMinorTick = 'on';
                    
                    
                    set(gca,'FontName','Times New Roman','FontSize',14,'LineWidth',2.5,'FontWeight','bold');
                    Tight = get(gca, 'TightInset');  %Gives you the bording spacing between plot box and any axis labels
                    %[Left Bottom Right Top] spacing
                    NewPos = [Tight(1) Tight(2) 1-Tight(1)-Tight(3)-0.01 1-Tight(2)-Tight(4)-0.01]; %New plot position [X Y W H]
                    set(gca, 'Position', NewPos);
                    %saveas(gcf,['/shared/dondo/home/ml423/world_mts_runs/figs/IRC_images/Basin' num2str(gid) '_' rainevent '_' num2str(ICn) '_W' num2str(itranum) num2str(wnum) '.jpg'])
                    %exportgraphics(gcf,['/taiga/ML/AI_20/pics/Basin' num2str(gid) '_ind_' num2str(bid) 'event' num2str(eid,'%2.2d') '_AI_20ens_' AItool '.jpg'],'Resolution',300)
                    %}
                end
                
                %}  
            end
            % done plot optimum Wxy
            
            
        end
    end
end
%%
andes_AIstats_tmp = andes_AIstats_RF;
ubid_andes = unique(andes_AIstats_tmp(:,1));
ubidinfo_andes = [];ubidfull_andes = {};
for i = 1:length(ubid_andes)
    tmp = find(andes_AIstats_tmp(:,1)==ubid_andes(i));
    AIori = nanmedian(andes_AIstats_tmp(tmp,3));
    AIopt = nanmedian(andes_AIstats_tmp(tmp,4));
    goodratio = length(find(andes_AIstats_tmp(tmp,4)>andes_AIstats_tmp(tmp,3)))/length(tmp);
    ubidinfo_andes(i,:) = [ubid_andes(i),AIori,AIopt,goodratio];
    ubidfull_andes{i} = [ones(length(tmp),1)*ubid_andes(i),andes_AIstats_tmp(tmp,3),andes_AIstats_tmp(tmp,4)];
    figure
    cdfplot(ubidfull_andes{i}(:,2))
    hold on
    cdfplot(ubidfull_andes{i}(:,3))
    title(num2str(ubidfull_andes{i}(1)))
end

locc = find(andes_AIstats_NN(1:3052,4)~=0&andes_AIstats_RF(1:3052,4)~=0);
figure
t = tiledlayout(1,1,'Padding','none');
t.Units = 'inches';
t.OuterPosition = [0.25 0.25 4.95 3.6];
nexttile;
h1 = cdfplot(andes_AIstats_NN(locc,3));
hold on
h2 = cdfplot(andes_AIstats_NN(locc,4));
hold on
h3 = cdfplot(andes_AIstats_RF(locc,4));
set(h1,'Color',[256 0 0]/256,'LineWidth',2)
set(h2,'Color',[75 68 56]/256,'LineWidth',2)
set(h3,'Color',[0 0 256]/256,'LineWidth',2)
ylabel('')
xlabel('')
title('AI performance','FontWeight','normal')
xlim([-3 1])
% hax = gca;
% hax.XAxis.MinorTickValues = linspace(0,200,41);
% hax.XMinorTick = 'on';
% xticks([0:10:100])
% xticklabels([0:10:100])
% hax.YAxis.MinorTickValues = linspace(0,1,11);
% hax.YMinorTick = 'on';
legend({'ERA5L','ERA5L-NN','ERA5L-RF'},'Location','Northwest')
set(gca,'FontName','Times New Roman','FontSize',16)%,'LineWidth',2.5,'FontWeight','bold');
Tight = get(gca, 'TightInset');  %Gives you the bording spacing between plot box and any axis labels
%[Left Bottom Right Top] spacing
NewPos = [Tight(1) Tight(2) 1-Tight(1)-Tight(3) 1-Tight(2)-Tight(4)]; %New plot position [X Y W H]
set(gca, 'Position', NewPos);% when use this set the figure visible off for better cropping figure

exportgraphics(gcf,['/shared/dondo/home/ml423/world_mts_runs/figs/synthesis/cross-regional-non-alps.jpg'],'Resolution',800)

%% plot AI outputs for Himalayas

clear tmpp
clear tmp
clear predatatmp
clear data3D
addpath('/shared/dondo/home/ml423/HMIOP/NoDAre/New era/')
addpath('/shared/dondo/home/ml423/DA/Critical matrices/')
addpath('/shared/dondo/home/ml423/DA/slanCM/')

tmp = load(['/shared/dondo/home/ml423/DA/Critical matrices/world_selected_gauge_events_allinfo_v3.mat'],'gauge_events_first_last');
Basins_sum = tmp.gauge_events_first_last;clear tmp;
tmp = load(['/shared/dondo/home/ml423/DA/Critical matrices/WORLD_events_v3.mat'],'WORLD_events');
WORLD_events = tmp.WORLD_events;clear tmp;

flowtype = {'streamflowarea','interflowarea','overlandflowarea','baseflowarea'};
ftp = 1;
AItool = 'A15';

aiens = 5;
hima_AIstats = [];

himaall = find(gauge_events_first_last(:,5)<1200);
hiall = [];
noflood = 442;
for ki = 1:length(himaall)
    bid = find(WORLD_events(:,1)==gauge_events_first_last(himaall(ki),1));
    if bid>=408&bid<=465
        if bid~=noflood
            if gauge_events_first_last(himaall(ki),3)>19900000
                hiall = [hiall;bid];
            end
        end
    end
end

hima_alleve = [];
for i = 1:length(hiall)
    bid = hiall(i);
    ixx = find(~isnan(WORLD_events(bid,:)));
    etot = length(ixx)-1;
    if etot == 0
        continue
    else
        if etot>149
            etot = 149
        end
        for eid = 1:etot
            close all
            % the if comments below are data quality control of the Alps
            % and hima
            %if isnan(event_quality(bid,eid+1))||isnan(kge_optm(bid,eid))
            if eqcontrol(bid,eid)==1%isempty(dssumtmp{bid,eid})||cvtmp(bid,eid)>1000
                hima_alleve = [hima_alleve;[bid,eid]];
            end
        end
    end
end


%andes_AIstats = [];
clear unfinish_hima
for ii = 1:length(hima_alleve)%[329 332 333 334 335 336 337 338 184 185 186 187 290 296 298 299 301 302 215 216 217 224 229 241 242 250 266 299 320 321 322 314]%[215 184]%[3 50]
    
    bid = hima_alleve(ii,1);
    gid = WORLD_events(bid,1);
    
    gloc = find(Basins_sum(:,1)==gid);
    
    r = Basins_sum(gloc,6);
    c = Basins_sum(gloc,7);
    rg = Basins_sum(gloc,8);
    cg = Basins_sum(gloc,9);
    area = Basins_sum(gloc,5);
    
    % if bid>=412
    %     cvtmp = cvfloodhima;
    %     dssumtmp = dssumhima;
    % elseif bid>=134
    %     cvtmp = cvfloodandes;
    %     dssumtmp = dssumandes;
    % end
    
    tmp = load(['/shared/dondo/home/ml423/DA/Critical matrices/World_Basins/inttind/Basin' num2str(gid) '_ind.mat'],'ind');
    ind = tmp.ind;
    tmp = load(['/shared/dondo/home/ml423/DA/Critical matrices/World_Basins/inttind/Basin' num2str(gid) '_intt.mat'],'intt');
    intt = tmp.intt;
    % ffnm=['/shared/dondo/home/ml423/world_mts/Basin' num2str(gid) '_fixed_900m/dem.bin'];
    % fid=fopen(ffnm,'rb','ieee-le');
    % tmpdata=fread(fid,inf,'single');
    % fclose(fid);
    % dem = reshape(tmpdata,c,r);
    % dem(ind) = nan;
    % demf = dem';
    % colo = slanCM('terrain',80);
    % colo(1,:) = [0.8 0.8 0.8];
    % figure
    % imagesc(demf)
    % colormap(colo)
    
    % ffnm=['/shared/dondo/home/ml423/world_mts/Basin' num2str(gid) '_fixed_900m/countfile.out'];
    % fid=fopen(ffnm);
    % facc=fscanf(fid,'%d',[c,r]);
    % fclose(fid);
    % facc(ind) = nan;
    % [streamx,streamy] = find(facc>5);
    
    
    
    ICn = 4140;
    ICC = 0;% IF ICC is part of the algo
    
    fileind = -100000*ICn;
    KGEsu = [];
    evesu = [];
    budgetalldbkc = [];
    budgetallirc = [];
    
    ixx = find(~isnan(WORLD_events(bid,:)));
    etot = length(ixx)-1;
    if etot == 0
        continue
    else
        if etot>149
           etot = 149 
        end
        for eidd = 1
            eid = hima_alleve(ii,2);
            %close all
            % the if comments below are data quality control of the Alps
            % and hima
            %if isnan(event_quality(bid,eid+1))||isnan(kge_optm(bid,eid))
            if eqcontrol(bid,eid)~=1%isempty(dssumtmp{bid,eid})||cvtmp(bid,eid)>1000    
                disp('Event quality is bad, peak too early or too late')
                continue
            else
                rainevent = num2str(WORLD_events(bid,eid+1));
                yth = floor(WORLD_events(bid,eid+1)/10000);
                dth = mod(WORLD_events(bid,eid+1),100);
                mth = floor((WORLD_events(bid,eid+1)-yth*10000)/100);
                
                thresdate = datetime(yth,mth,dth);
                
                durati = (thresdate-days(15)):days(1):(thresdate+days(15));
                [dury,durm,durd] = ymd(durati);
                tst1 = num2str(durm(1)*100+durd(1),'%4.4d');
                tst2 = num2str(durm(30)*100+durd(30),'%4.4d');
                
                xtlabels ={};
                for tl = 1:30
                    xtlabels{tl} =  [ num2str(durm(tl),'%2.2d') '/' num2str(durd(tl),'%2.2d') ];
                end
                
                
                
                
                
                ffnm = ['/shared/dondo/home/ml423/world_mts_runs/Basins_events_chars/'];
                fid = fopen([ffnm 'Basin' num2str(gid) '_event' rainevent '.txt'],'r');
                if fid == -1
                    %fclose(fid);
                    continue
                end
                windowtimes = fscanf(fid,'%d\n',[4,1]);
                fclose(fid);
                if length(windowtimes)<4
                   continue 
                end
                %close all
                
                
                if 1==1
                    
                    %fh = APL_frontcut(bid,eid);
                    
                    tstr = [num2str(dury(1)) '-' num2str(durm(1),'%2.2d') '/' num2str(durd(1),'%2.2d') '-' num2str(durm(30),'%2.2d') '/' num2str(durd(30),'%2.2d')];
                    
                    
                    obstep = 90;
                    
                    tmp = load(['/shared/dondo/home/ml423/world_mts/obs_events/events/obs_' num2str(gid) '_event' rainevent '.mat']);
                    tmp = tmp.obsnow;
                    obsnow = tmp;
                    obsnan=find(isnan(obsnow));
                    if isempty(obsnan)
                    else
                        for na = 1:length(obsnan)
                            obsnow(obsnan(na)) = 0.5*(obsnow(obsnan(na)-1)+obsnow(obsnan(na)+1));
                        end
                    end
                    obs2 = interp1(1:31,obsnow,[1:0.33333:31]);
                    obsnow = obs2(1:90);
                    
                    %figure
                    %plot(obsnow)
                    
                    tmp = load(['/shared/dondo/home/ml423/world_mts/Precip/Event/ERA5_Basin' num2str(gid) '_' rainevent '.mat']);
                    s4dbkc = tmp.rain;
                    
                    MCp = zeros(c,r);
                    pts = [];
                    dbkcr = 0;rmax = [];
                    for i = 1:720
                        % MCp = MCp + s4dbkc{i};
                        % pts(i) = mean(s4dbkc{i}(intt));
                        % dbkcr = dbkcr+pts(i);
                        rmax(i) = max(s4dbkc{i}(intt));
                    end
                    % MCp(ind) = nan;
                    % ori_mean = nanmean(MCp(:));
                    % ama = jet(40);ama(1,:) = [1 1 1];
                    % figure
                    % imagesc(MCp')
                    % colorbar
                    % colormap(ama)
                    % caxis([0.2*min(MCp(:)) 1.1*max(MCp(:))])
                    % 
                    
                    
                    % for Nature Science paper, consider the major events only
                    
                    launchp = find(rmax>0.5);% when rainfall >1mm/hr launch particles
                    ptmaxthre = 0.5*max(rmax);% used as a threshold to separate events
                    la2 = find(launchp<(windowtimes(4)*8));
                    la1 = find(launchp<(windowtimes(3)*8));
                    accu0 = 0;
                    if isempty(la1)
                        continue
                    end
                    for tmpla = launchp(la1(end)):-1:1
                        if ismember(tmpla,launchp)
                            accu0 = 0;
                            rainbegintime = tmpla;
                        else
                            accu0 = accu0+1;
                        end
                        
                        if area>1000
                            noraininterval = 12;% for big basins, it is likely there is always rain in the basin, therefore using a small norain limit to distinguish different rain events.
                        else
                            noraininterval = 24;
                        end
                        if accu0>noraininterval
                            % consecutive x hours of basin max rain less than 0.5mm/hr,
                            % this should be counted as a different rain cell, should not
                            % be tracked!
                            % one exception is if another major rain is 48-24 ahead, should be
                            % considered as well
                            rsmallest = max([(tmpla-noraininterval),1]);
                            
                            reconsider_period = rsmallest:tmpla;
                            %reconsider_period = (tmpla-noraininterval):tmpla;
                            secondary_rain = find(rmax(reconsider_period)>0.5);
                            
                            if max(rmax(reconsider_period))>ptmaxthre
                                % adjust the rainbegintime here
                                rainbegintime = reconsider_period(secondary_rain(1));
                            end
                            smalltmpla = max([floor(tmpla/8),2]);
                            if obsnow(smalltmpla)>obsnow(smalltmpla-1)
                                % meaning obsflow is still rising.
                                if ~isempty(secondary_rain)
                                    rainbegintime = reconsider_period(secondary_rain(1));
                                end
                            end
                            
                            break
                        end
                    end
                    
                    la0 = find(launchp == rainbegintime);
                    
                    refinedlaunchp = launchp(la0:la2(end));
                    stat_left_bd_hydro = windowtimes(1);% systematic defined
                    stat_left_bd_rain = floor(rainbegintime/8);%target rainevent adjusted
                    stat_left_bd = max([stat_left_bd_hydro,stat_left_bd_rain]);
                    
                    ICtiming = refinedlaunchp(1);%1-720, can be used to extract original IC from IRCICC studies
                    
                    upcap = min([90 (windowtimes(4)+6)]);
                    
                    if upcap<=stat_left_bd
                        continue
                    end
                    
                    no3sum = nan(aiens,720);
                    %no2sum = [];
                    rts1 = eval_ind{bid,eid};
                    AI_mean_rain = zeros(c,r,length(rts1));
                    for replic = 1:aiens
                        
                        ffnmx=['/taiga/ML/AI_20/Basin' num2str(gid) 'outputs/' rainevent '_AI_' num2str(replic,'%2.2d') '_' AItool '/' rainevent '/' flowtype{ftp}];
                        fid=fopen(ffnmx,'rb','ieee-le');
                        if fid == -1
                            disp('no AI done')
                            aibreak = 1;
                            break
                        else
                            aibreak = 0;
                        end
                        tmpdata=fread(fid,inf,'single');
                        fclose(fid);
                        if length(tmpdata)<(c*r*720)
                            disp([ii,ii])
                            unfinish_hima(ii,:) = [bid,eid];
                            aibreak = 1;
                            break
                        end
                        NODAc2o5=reshape(tmpdata,c,r,720);
                        no3 = [];
                        for k = 1:720
                            no3(k) = abs(NODAc2o5(cg,rg,k));
                        end
                       
                        no3sum(replic,:) = no3;
                        % ffnm=['/taiga/ML/AI_20/Basin' num2str(gid) 'outputs/' rainevent '_AI_' num2str(replic,'%2.2d') '_' AItool '/' rainevent '/precipitation'];
                        % fid=fopen(ffnm,'rb','ieee-le');
                        % tmpdata=fread(fid,inf,'single');
                        % fclose(fid);
                        % NODAc2o5=reshape(tmpdata,c,r,720);
                        % flownoda5 = [];
                        % meanr250 = 0;
                        % tmpsum = zeros(c,r);
                        % for k = 1:720
                        %     tmp = NODAc2o5(:,:,k);
                        %     tmpsum = tmpsum+tmp;
                        %     flownoda5(k) = mean(tmp(intt));
                        %     meanr250 = meanr250+flownoda5(k);
                        % end
                        % no2 = flownoda5;%GREEN
                        % rainsd = std(tmpsum(intt));
                        % no2sum = [no2sum;no2];
                        % 
                        % if replic == 1
                        %     figure
                        %     imagesc(tmpsum')
                        %     colorbar
                        %     colormap(jet(30))
                        % end
                        ffnmx=['/shared/dondo/home/ml423/world_mts/Basin' num2str(gid) '_input_ERA5_900m_tmp/' rainevent '/Precipitation_ERA5_lRFR_ori_' num2str(replic,'%3.3d') ];
                        fid=fopen(ffnmx,'rb','ieee-le');
                        tmpdata=fread(fid,inf,'single');
                        fclose(fid);
                        NODAc2o5=reshape(tmpdata,c,r,720);
                        
                        AI_mean_rain = AI_mean_rain+NODAc2o5(:,:,rts1);

                    end
                    
                    if aibreak==1
                       continue 
                    end

                    rain_AI{bid,eid} = AI_mean_rain/aiens;
                    
                    no3vol_rank = sum(no3sum,2);
                    [~,mloc] = sort(no3vol_rank);
                    
                    
                    if aiens>5
                        no3m = mean(no3sum(mloc(6:15),:),1);   
                    else
                        no3m = mean(no3sum,1);
                    end
                    
%                     filenam = str2num(rainevent)*100+fileind;
%                     ffnmx=['/shared/dondo/home/ml423/world_mts_runs/IRC1/Basin' num2str(gid) 'outputs/refxtmp' num2str(filenam) 'fullresults/' rainevent '_SW/' rainevent '/' flowtype{ftp}];
%                     fid=fopen(ffnmx,'rb','ieee-le');
%                     tmpdata=fread(fid,inf,'single');
%                     fclose(fid);
%                     NODAc2o5=reshape(tmpdata,c,r,720);
%                     flownoda5 = [];
%                     for k = 1:720
%                         flownoda5(k) = abs(NODAc2o5(cg,rg,k));
%                     end
%                     no5 = flownoda5;% original era5 with original IC
%                     
                    ffnmx=['/taiga/ML/AI_20/Basin' num2str(gid) 'outputs/' rainevent '_AI_' num2str(0,'%2.2d') '/' rainevent '/' flowtype{ftp}];
                    fid=fopen(ffnmx,'rb','ieee-le');
                    if fid == -1
                        continue
                    end
                    tmpdata=fread(fid,inf,'single');
                    fclose(fid);
                    NODAc2o5=reshape(tmpdata,c,r,720);
                    flownoda5 = [];
                    for k = 1:720
                        flownoda5(k) = abs(NODAc2o5(cg,rg,k));
                    end
                    no6 = flownoda5;% original era5 with best IC from IRCicc
                    
                    
                    cal_window = (stat_left_bd-1)*8+1:(upcap-1)*8+1;
                    % mass balance check fo W00 and W34
                    %{
                        filenamp = str2num(rainevent)*100+fileind+2;% the very first P
                        cc = ls(['/shared/dondo/home/ml423/HMIOP/NoDAre/IRC_LB/Basin' num2str(bid,'%2.2d') 'outputs/refxtmp' num2str(filenamp) '/xtmpbb1*']);
                        cc2 = strsplit(cc);
                        tmpC  = struct2cell(load([cc2{end-1}]));
                        P00 = tmpC{1};
                        filenamp = str2num(rainevent)*100+fileind+24;% the very first P
                        cc = ls(['/shared/dondo/home/ml423/HMIOP/NoDAre/IRC_LB/Basin' num2str(bid,'%2.2d') 'outputs/refxtmp' num2str(filenamp) '/xtmpbb1*']);
                        cc2 = strsplit(cc);
                        tmpC  = struct2cell(load([cc2{1}]));
                        P34 = tmpC{1};
                        P00tot = 0;
                        P34tot = 0;
                        for k = 1:288
                            tmp = P00{k};
                            P00tot = P00tot+mean(tmp(intt));
                            tmp = P34{k};
                            P34tot = P34tot+mean(tmp(intt));
                        end
                        % mass balance quick calculations
                        runoff250 = sum(no3)*300/length(intt)/250/250*1000;
                        
                        
                        
                        fixedstr = ['/shared/dondo/home/ml423/HMIOP/NoDAre/LB/APL_Basin' num2str(bid,'%2.2d') '_fixed_250m/'];
                        datastr = ['/shared/dondo/home/ml423/HMIOP/NoDAre/IRC_LB/Basin' num2str(bid,'%2.2d') 'outputs/refxtmp' num2str(filenam) 'fullresults/Basin' num2str(bid,'%2.2d') '_output_250m5min_' rainevent '_SW/' rainevent '00/'];
                        
                        lastnote = 'c';
                        if p==0
                            laststr = ['/shared/dondo/home/ml423/HMIOP/NoDAre/LB/Basin' num2str(bid,'%2.2d') '_output_250m5min_' rainevent '_' lastnote '/'];
                            [smdelta250,gddelta250] = waterbalance(0,fixedstr,c,r,288,288,laststr,datastr,'',intt,250);
                        elseif p==9
                            if ICC==1% if ICC is part of the algo
                                laststr = ['/shared/dondo/home/ml423/HMIOP/NoDAre/IRC_LB/ICall/' num2str(ICn) '/Basin' num2str(bid,'%2.2d') '_output_250m5min_' rainevent lastnote '/'];
                            else
                                laststr = ['/shared/dondo/home/ml423/HMIOP/NoDAre/LB/Basin' num2str(bid,'%2.2d') '_output_250m5min_' rainevent '_' lastnote '/'];
                            end
                            [smdelta250,gddelta250] = waterbalance(0,fixedstr,c,r,288,288,laststr,datastr,'',intt,250);
                        end
                        
                        if p==0
                            rainDBKC250 = dbkcr;% in mm
                            riverDBKC250 = runoff250;% in mm
                            deltasDBKC250 = smdelta250/250/250/length(intt)*1000; %in mm
                            gdDBKC250 = gddelta250/250/250/length(intt)*1000; %in mm
                            budgetDBKC250 = rainDBKC250-riverDBKC250-(deltasDBKC250+gdDBKC250);
                            wbterms_DBKC = [rainDBKC250,riverDBKC250,deltasDBKC250,gdDBKC250,budgetDBKC250];
                            budgetalldbkc = [budgetalldbkc;[str2num(rainevent),wbterms_DBKC]];
                        elseif p==9
                            rainIRC250 = P34tot;% in mm
                            riverIRC250 = runoff250;% in mm
                            deltasIRC250 = smdelta250/250/250/length(intt)*1000; %in mm
                            gdIRC250 = gddelta250/250/250/length(intt)*1000; %in mm
                            budgetIRC250 = rainIRC250-riverIRC250-(deltasIRC250+gdIRC250);
                            wbterms_IRC = [rainIRC250,riverIRC250,deltasIRC250,gdIRC250,budgetIRC250];
                            budgetallirc = [budgetallirc;[str2num(rainevent),wbterms_IRC]];
                        end
                        
                        
                    %}
                    
                    
                    simu = no6(1:8:end);
                    
                    
                    simustat = simu(stat_left_bd:upcap);
                    obsnowstat = obsnow(stat_left_bd:upcap);
                    tmp = 0;
                    tmp1 = 0;
                    tmp2 = mean(obsnowstat);
                    for i = 1:length(obsnowstat)
                        tmp = tmp + (simustat(i)-obsnowstat(i))^2;
                        tmp1 = tmp1 + (tmp2-obsnowstat(i))^2;
                    end
                    NSEori = 1-tmp./tmp1;
                    
                    rr = corrcoef(simustat',obsnowstat);
                    rstar = rr(1,2);
                    obsstd = std(obsnowstat);
                    simustd = std(simustat);
                    obsmean = mean(obsnowstat);
                    simumean = mean(simustat);
                    KGEori = 1-sqrt((rstar-1)^2+(simustd/obsstd-1)^2+(simumean/obsmean-1)^2);
                    
                    simu = no3m(1:8:end);
                    simustat = simu(stat_left_bd:upcap);
                    obsnowstat = obsnow(stat_left_bd:upcap);
                    tmp = 0;
                    tmp1 = 0;
                    tmp2 = mean(obsnowstat);
                    for i = 1:length(obsnowstat)
                        tmp = tmp + (simustat(i)-obsnowstat(i))^2;
                        tmp1 = tmp1 + (tmp2-obsnowstat(i))^2;
                    end
                    NSE = 1-tmp./tmp1;
                    
                    rr = corrcoef(simustat',obsnowstat);
                    rstar = rr(1,2);
                    obsstd = std(obsnowstat);
                    simustd = std(simustat);
                    obsmean = mean(obsnowstat);
                    simumean = mean(simustat);
                    KGE = 1-sqrt((rstar-1)^2+(simustd/obsstd-1)^2+(simumean/obsmean-1)^2);
              
                    [vsimu,tsimu]=max(simustat);
                    [vobs,tobs]=max(obsnowstat);
                    EPV = vsimu-vobs;
                    EPT = (tsimu-tobs)*1;
                    simvol = 0;
                    obsvol = 0;
                    for vlo=2:length(obsnowstat)-1
                        obsvol = obsvol+obsnowstat(vlo)*1;
                        simvol = simvol+simustat(vlo)*1;
                    end
                    obsvol = obsvol+0.5*(obsnowstat(1)+obsnowstat(end))*1;
                    simvol = simvol+0.5*(simustat(1)+simustat(end))*1;
                    
                    EV = (simvol-obsvol)./obsvol;
                    hima_AIstats(ii,:) = [bid,eid,KGEori,KGE];
                    %{
                    % KGE is period defined, while ther other stats
                    % use the whole series
                    
                    %                             disp(['NSE is ' num2str(NSE)])
                    %                             disp(['KGE is ' num2str(KGE)])
                    %                             disp(['EPV is ' num2str(EPV)])
                    %                             disp(['EPT is ' num2str(EPT)])
                    
                    
%                     
%                                                 figure
%                                                 t = tiledlayout(1,1,'Padding','none');
%                                                 t.Units = 'inches';
%                                                 t.OuterPosition = [0.25 0.25 4.8 3.2];
%                                                 nexttile;
                    
                    f = figure('visible','on');
                    f.Position = [5 5 700 350];
%                     
                    pos1 = [0.13 0.8 0.67 0.15];
                    %subplot('Position',pos1)
                    pos2 = [0.15 0.32 0.78 0.42];
                    
                    subplot('Position',pos2)
                    % hjh = plot(2:3:obstep*3,obsdatav1(:,6),'k',1:288,no3,'g',1:288,no2,'r',1:288,no5,'m');
                    % set(hjh,'LineWidth',2.5)
                    % legend('Obs.','IMERG_D','STIV_D','STIV_{DBK}','Location','Northeast')
                    %yyaxis left
                    xj = 2:8:obstep*8;
                    xj = xj';
                    bd1 = 1.05*obsnow;
                    bd2 = 0.95*obsnow;
                    sh = scatter(2:8:obstep*8,obsnow,45,'g','filled');
                    sh.MarkerEdgeColor = [0,1,0];
                    for ii = 1:size(no3sum,1)
                        hold on
                        hjh = plot(1:8:720,no3sum(ii,1:8:720),'-.','LineWidth',1.5,'Color',[0, 0, 1, 0.3]);
                    end
                   
%                     hold on
%                     hjh = plot(1:8:720,no5(1:8:720),'m:','LineWidth',2);
                    hold on
                    hjh = plot(1:8:720,no6(1:8:720),'r:','LineWidth',2);
                    hold on
                    hjh = plot(1:8:720,no3m(1:8:720),'b-.');
                    hold on
                    hjb = plot([1+8*stat_left_bd 1+8*stat_left_bd],[0 1000000],'k-.');
                    hold on
                    hjb2 = plot([1+8*(windowtimes(4)+6) 1+8*(windowtimes(4)+6)],[0 1000000],'k-.');
                    
                    % ymax is based on max of y but rounded to
                    % nearest interval defined by max of y
                   
                        ymm = ceil(max(max(1.5*obsnow),max(1.5*no3m(:)))/(max(obsnow)/4))*(max(obsnow)/4);
                    
                    tf = text(16,0.85*ymm,['' num2str(KGEori,'%2.2f')],'Color','r');
                    tf.FontSize = 14;
                    tf.FontWeight = 'bold';
                    
                    tf = text(16,0.72*ymm,['' num2str(KGE,'%2.2f')],'Color','b');
                    tf.FontSize = 14;
                    tf.FontWeight = 'bold';
                    set(hjh,'LineWidth',2.5)
                    %             hold on
                    %             patch([xj; flipud(xj)]', [bd1; flipud(bd2)]',  'k', 'FaceAlpha',0.2, 'EdgeColor','none')
                    %             hold off
                    ylim([0 ymm])
                    ylabel('Streamflow (m^3/s)')
                    
                    box on
                    
                    %title([num2str(gid) ' IRC W' num2str(itranum) num2str(wnum) ' | ' num2str(area,'%0.0f') ' km^2'])
                    
                    
                    xlabel([tstr ' '])
                    %set(gca,'XTick',1:3:288,'XTickLabel',[])
                    xlim([0 719])
                    set(gca,'XTick',0:144:720,'XTickLabel',xtlabels(1:6:30),'FontSize',12,'LineWidth',2.5,'FontWeight','bold')
                    
                    
                    hax = gca;
                    hax.XAxis.MinorTickValues = linspace(0,720,31);
                    hax.XMinorTick = 'on';
                    
                    
                    set(gca,'FontName','Times New Roman','FontSize',14,'LineWidth',2.5,'FontWeight','bold');
                    Tight = get(gca, 'TightInset');  %Gives you the bording spacing between plot box and any axis labels
                    %[Left Bottom Right Top] spacing
                    NewPos = [Tight(1) Tight(2) 1-Tight(1)-Tight(3)-0.01 1-Tight(2)-Tight(4)-0.01]; %New plot position [X Y W H]
                    set(gca, 'Position', NewPos);
                    %saveas(gcf,['/shared/dondo/home/ml423/world_mts_runs/figs/IRC_images/Basin' num2str(gid) '_' rainevent '_' num2str(ICn) '_W' num2str(itranum) num2str(wnum) '.jpg'])
                    %exportgraphics(gcf,['/taiga/ML/AI_20/pics/Basin' num2str(gid) '_ind_' num2str(bid) 'event' num2str(eid,'%2.2d') '_AI_20ens_' AItool '.jpg'],'Resolution',300)
                    %}
                end
                
                %}  
            end
            % done plot optimum Wxy
            
            
        end
    end
end

ubid_hima = unique(hima_AIstats(:,1));
ubidinfo_hima = [];ubidfull_hima = {};
for i = 1:length(ubid_hima)
    tmp = find(hima_AIstats(:,1)==ubid_hima(i));
    AIori = nanmedian(hima_AIstats(tmp,3));
    AIopt = nanmedian(hima_AIstats(tmp,4));
    goodratio = length(find(hima_AIstats(tmp,4)>hima_AIstats(tmp,3)))/length(tmp);
    ubidinfo_hima(i,:) = [ubid_hima(i),AIori,AIopt,goodratio];
    ubidfull_hima{i} = [ones(length(tmp),1)*ubid_hima(i),hima_AIstats(tmp,3),hima_AIstats(tmp,4)];
    figure
    cdfplot(ubidfull_hima{i}(:,2))
    hold on
    cdfplot(ubidfull_hima{i}(:,3))
    title(num2str(ubidfull_hima{i}(1)))
end
figure
cdfplot(hima_AIstats(:,3))
hold on
cdfplot(hima_AIstats(:,4))
%% [key plots] plot Alps, Andes, Himalaya AI results
% use the readhgt.m file in /DA/Critical matrices/hgt_dem_data/ folder so
% that it is gray scale.
addpath('/shared/dondo/home/ml423/DA/Critical matrices/')

Mb = shaperead('world-administrative-boundaries');
kge_ai_ratio = [];
bai = unique(alps_AIstats(:,1));

for i = 1:length(bai)
    if bai(i)==0
        continue
    end
    
    orik = alps_AIstats(find(alps_AIstats(:,1)==bai(i)),3);
    optk = alps_AIstats(find(alps_AIstats(:,1)==bai(i)),4);
    kratio = length(find(optk>orik))/length(optk);
    if kratio>=0.8
        kgesym = 0.9;
    elseif kratio>=0.6
        kgesym = 0.7;
    elseif kratio>=0.4
        kgesym = 0.5;
    elseif kratio>=0.2
        kgesym = 0.3;
    else
        kgesym = 0.1;
    end
    kge_ai_ratio(i,:) = [outletlatlon(bai(i),1:3),kgesym];
end

AA = jet(30);
AB = [AA(2,:);AA(10,:);[1 1 1];AA(23,:);AA(29,:)];

% lat lon range of the Alps
alpslat1 = 43.5;
alpslat2 = 49;
alpslon1 = 4;
alpslon2 = 17;

% ori ERA5 amounts

figure
t = tiledlayout(1,1,'Padding','none');
t.Units = 'inches';
t.OuterPosition = [0.25 0.25 5.1 4.2];
nexttile;
readhgt(alpslat1:alpslat2,alpslon1:alpslon2)
alpha(.2)
title([' '],'FontWeight','normal')
hold on 
for k=1:length(Mb)
     line(Mb(k).X(:),Mb(k).Y(:),'Color','k','LineWidth',0.5,'LineStyle',':'); 
end
hold on
sc = scatter(kge_ai_ratio(:,3),kge_ai_ratio(:,2),26,kge_ai_ratio(:,4),'filled');
sc.MarkerEdgeColor = 'k';
sc.LineWidth = 1;
ylim([alpslat1 alpslat2])
xlim([alpslon1 alpslon2])
colormap([AB])
clim([0 1])
cb = colorbar;
a = cb.Position;
set(cb,'Position',[a(1)+0.03 a(2)-0.06 0.02 0.65])
set(gca,'FontName','Times New Roman','FontSize',14)%,'LineWidth',2.5);%,'FontWeight','bold');
box on  
Tight = get(gca, 'TightInset');  %Gives you the bording spacing between plot box and any axis labels
%[Left Bottom Right Top] spacing
NewPos = [Tight(1) Tight(2) 1-Tight(1)-Tight(3)-0.15 1-Tight(2)-Tight(4)]; %New plot position [X Y W H]
set(gca, 'Position', NewPos);
                          
exportgraphics(gcf,['/shared/dondo/home/ml423/world_mts_runs/figs/synthesis/alps_AIRF_stats.jpg'],'Resolution',500)



figure
t = tiledlayout(1,1,'Padding','none');
t.Units = 'inches';
t.OuterPosition = [0.25 0.25 4.95 3.6];
nexttile;
cd1 = cdfplot(alps_AIstats(alps_AIstats(:,1)~=0,3));set(cd1,'Color','r','LineWidth',2);
hold on
cd2 = cdfplot(alps_AIstats(alps_AIstats(:,1)~=0,4));set(cd2,'Color','b','LineWidth',2);
hold off
xlim([-3 1])
ylim([0 1.02])
xlabel('KGE')
ylabel('CDF')
title(' ')
box on
legend(["ERA5L" "ERA5L-AI"],'Location','Northwest')
set(gca,'FontName','Times New Roman','FontSize',14)
exportgraphics(gcf,['/shared/dondo/home/ml423/world_mts_runs/figs/synthesis/KGE_CDF_comparison_AIRF.jpg'],'Resolution',600)

%% Andes
kge_ai_ratio_andes = [];
bai = unique(andes_AIstats(:,1));
for i = 1:length(bai)
    if bai(i)==0
        continue
    end
    
    orik = andes_AIstats(find(andes_AIstats(:,1)==bai(i)),3);
    optk = andes_AIstats(find(andes_AIstats(:,1)==bai(i)),4);
    kratio = length(find(optk>orik))/length(optk);
    if kratio>=0.8
        kgesym = 0.9;
    elseif kratio>=0.6
        kgesym = 0.7;
    elseif kratio>=0.4
        kgesym = 0.5;
    elseif kratio>=0.2
        kgesym = 0.3;
    else
        kgesym = 0.1;
    end
    kge_ai_ratio_andes(i,:) = [outletlatlon(bai(i),1:3),kgesym,bai(i)];
end
AA = jet(30);
AB = [AA(2,:);AA(10,:);[1 1 1];AA(23,:);AA(29,:)];

% lat lon range of the Andes
andes2lat1 = -6;
andes2lat2 = 24;
andes2lon1 = -95;
andes2lon2 = -60;
%
areamap=1;
cdfmap=1;
for andeszoom = 1:4
    if andeszoom == 1
        %entire South America
        andes2lat1 = -43;
        andes2lat2 = 24;
        andes2lon1 = -95;
        andes2lon2 = -34;
        nts = 'SA';
        basinloc = find(kge_ai_ratio_andes(:,2)>andes2lat1&kge_ai_ratio_andes(:,2)<andes2lat2&kge_ai_ratio_andes(:,3)>andes2lon1&kge_ai_ratio_andes(:,3)<andes2lon2);
        basinlist = kge_ai_ratio_andes(basinloc,5);
    elseif andeszoom == 2
        %pahamas
        andes2lat1 = 7;
        andes2lat2 = 13.8;
        andes2lon1 = -88;
        andes2lon2 = -78;
        nts = 'PAHAMA';
        basinloc = find(kge_ai_ratio_andes(:,2)>andes2lat1&kge_ai_ratio_andes(:,2)<andes2lat2&kge_ai_ratio_andes(:,3)>andes2lon1&kge_ai_ratio_andes(:,3)<andes2lon2);
        basinlist = kge_ai_ratio_andes(basinloc,5);

    elseif andeszoom == 3
        % Jamaca
        andes2lat1 = 17.5;
        andes2lat2 = 18.8;
        andes2lon1 = -78.5;
        andes2lon2 = -75.6;
        nts = 'JAMACA';
        basinloc = find(kge_ai_ratio_andes(:,2)>andes2lat1&kge_ai_ratio_andes(:,2)<andes2lat2&kge_ai_ratio_andes(:,3)>andes2lon1&kge_ai_ratio_andes(:,3)<andes2lon2);
        basinlist = kge_ai_ratio_andes(basinloc,5);

    elseif andeszoom == 4
        %puertorico
        andes2lat1 = 17.7;
        andes2lat2 = 18.7;
        andes2lon1 = -67.5;
        andes2lon2 = -65.4;
        nts = 'PRICO';
        basinloc = find(kge_ai_ratio_andes(:,2)>andes2lat1&kge_ai_ratio_andes(:,2)<andes2lat2&kge_ai_ratio_andes(:,3)>andes2lon1&kge_ai_ratio_andes(:,3)<andes2lon2);
        basinlist = kge_ai_ratio_andes(basinloc,5);
    end


    if areamap==1
        figure
        t = tiledlayout(1,1,'Padding','none');
        t.Units = 'inches';
        t.OuterPosition = [0.25 0.25 5.1 4.2];
        nexttile;
        readhgt(andes2lat1:andes2lat2,andes2lon1:andes2lon2,'srtm3')
        alpha(.2)
        title([' '],'FontWeight','normal')
        hold on
        for k=1:length(Mb)
            line(Mb(k).X(:),Mb(k).Y(:),'Color','k','LineWidth',0.5,'LineStyle',':');
        end
        hold on
        sc = scatter(kge_ai_ratio_andes(:,3),kge_ai_ratio_andes(:,2),26,kge_ai_ratio_andes(:,4),'filled');
        sc.MarkerEdgeColor = 'k';
        sc.LineWidth = 1;
        ylim([andes2lat1 andes2lat2])
        xlim([andes2lon1 andes2lon2])
        colormap([AB])
        clim([0 1])
        cb = colorbar;
        a = cb.Position;
        set(cb,'Position',[a(1)+0.03 a(2)-0.06 0.02 0.65])
        set(gca,'FontName','Times New Roman','FontSize',14)%,'LineWidth',2.5);%,'FontWeight','bold');
        box on
        Tight = get(gca, 'TightInset');  %Gives you the bording spacing between plot box and any axis labels
        %[Left Bottom Right Top] spacing
        NewPos = [Tight(1) Tight(2) 1-Tight(1)-Tight(3)-0.15 1-Tight(2)-Tight(4)]; %New plot position [X Y W H]
        set(gca, 'Position', NewPos);
        exportgraphics(gcf,['/shared/dondo/home/ml423/world_mts_runs/figs/synthesis/andes_AIRF_stats_' nts '.jpg'],'Resolution',700)
    end

    if cdfmap == 1
        locc = ismember(andes_AIstats(:,1),basinlist);
        figure
        t = tiledlayout(1,1,'Padding','none');
        t.Units = 'inches';
        t.OuterPosition = [0.25 0.25 4.95 3.6];
        nexttile;
        cd1 = cdfplot(andes_AIstats(locc==1,3));set(cd1,'Color','r','LineWidth',2);
        hold on
        cd2 = cdfplot(andes_AIstats(locc==1,4));set(cd2,'Color','b','LineWidth',2);
        hold off
        xlim([-3 1])
        ylim([0 1.02])
        xlabel('KGE')
        ylabel('CDF')
        title(' ')
        box on
        %legend(["ERA5L" "ERA5L-AI"],'Location','Northwest')
        set(gca,'FontName','Times New Roman','FontSize',14)
        if andeszoom>1
            set(gca,'FontName','Times New Roman','FontSize',16)
        end
        exportgraphics(gcf,['/shared/dondo/home/ml423/world_mts_runs/figs/synthesis/KGE_CDF_comparison_AIRF_Andes_' nts '.jpg'],'Resolution',600)
    end

end
        

%% Himalayas
kge_ai_ratio_hima = [];
bai = unique(hima_AIstats(:,1));
for i = 1:length(bai)
    if bai(i)==0
        continue
    end
    
    orik = hima_AIstats(find(hima_AIstats(:,1)==bai(i)),3);
    optk = hima_AIstats(find(hima_AIstats(:,1)==bai(i)),4);
    kratio = length(find(optk>orik))/length(optk);
    if kratio>=0.8
        kgesym = 0.9;
    elseif kratio>=0.6
        kgesym = 0.7;
    elseif kratio>=0.4
        kgesym = 0.5;
    elseif kratio>=0.2
        kgesym = 0.3;
    else
        kgesym = 0.1;
    end
    kge_ai_ratio_hima(i,:) = [outletlatlon(bai(i),1:3),kgesym,bai(i)];
end

areamap=1;
cdfmap=1;
for himazoom = 1:2
    if himazoom == 1
        % Southeast asia
        hima2lat1 = 15;
        hima2lat2 = 22;
        hima2lon1 = 98;
        hima2lon2 = 107;
        nts = 'Asia';
        basinloc = find(kge_ai_ratio_hima(:,2)>hima2lat1&kge_ai_ratio_hima(:,2)<hima2lat2&kge_ai_ratio_hima(:,3)>hima2lon1&kge_ai_ratio_hima(:,3)<hima2lon2);
        basinlist = kge_ai_ratio_hima(basinloc,5);
    elseif himazoom == 2
        %pahamas
        hima2lat1 = 26;
        hima2lat2 = 30;
        hima2lon1 = 82;
        hima2lon2 = 88;
        nts = 'Hima';
        basinloc = find(kge_ai_ratio_hima(:,2)>hima2lat1&kge_ai_ratio_hima(:,2)<hima2lat2&kge_ai_ratio_hima(:,3)>hima2lon1&kge_ai_ratio_hima(:,3)<hima2lon2);
        basinlist = kge_ai_ratio_hima(basinloc,5);

    elseif himazoom == 3
        % Jamaca
        hima2lat1 = 17.5;
        hima2lat2 = 18.8;
        hima2lon1 = -78.5;
        hima2lon2 = -75.6;
        nts = 'JAMACA';
        basinloc = find(kge_ai_ratio_hima(:,2)>hima2lat1&kge_ai_ratio_hima(:,2)<hima2lat2&kge_ai_ratio_hima(:,3)>hima2lon1&kge_ai_ratio_hima(:,3)<hima2lon2);
        basinlist = kge_ai_ratio_hima(basinloc,5);

    elseif himazoom == 4
        %puertorico
        hima2lat1 = 17.7;
        hima2lat2 = 18.7;
        hima2lon1 = -67.5;
        hima2lon2 = -65.4;
        nts = 'PRICO';
        basinloc = find(kge_ai_ratio_hima(:,2)>hima2lat1&kge_ai_ratio_hima(:,2)<hima2lat2&kge_ai_ratio_hima(:,3)>hima2lon1&kge_ai_ratio_hima(:,3)<hima2lon2);
        basinlist = kge_ai_ratio_hima(basinloc,5);
    end


    if areamap==1
        figure
        t = tiledlayout(1,1,'Padding','none');
        t.Units = 'inches';
        t.OuterPosition = [0.25 0.25 5.1 4.2];
        nexttile;
        readhgt(hima2lat1:hima2lat2,hima2lon1:hima2lon2,'srtm3')
        alpha(.2)
        title([' '],'FontWeight','normal')
        hold on
        for k=1:length(Mb)
            line(Mb(k).X(:),Mb(k).Y(:),'Color','k','LineWidth',0.5,'LineStyle',':');
        end
        hold on
        sc = scatter(kge_ai_ratio_hima(:,3),kge_ai_ratio_hima(:,2),46,kge_ai_ratio_hima(:,4),'filled');
        sc.MarkerEdgeColor = 'k';
        sc.LineWidth = 1;
        ylim([hima2lat1 hima2lat2])
        xlim([hima2lon1 hima2lon2])
        colormap([AB])
        clim([0 1])
        cb = colorbar;
        a = cb.Position;
        set(cb,'Position',[a(1)+0.03 a(2)-0.06 0.02 0.65])
        set(gca,'FontName','Times New Roman','FontSize',14)%,'LineWidth',2.5);%,'FontWeight','bold');
        box on
        Tight = get(gca, 'TightInset');  %Gives you the bording spacing between plot box and any axis labels
        %[Left Bottom Right Top] spacing
        NewPos = [Tight(1) Tight(2) 1-Tight(1)-Tight(3)-0.15 1-Tight(2)-Tight(4)]; %New plot position [X Y W H]
        set(gca, 'Position', NewPos);
        exportgraphics(gcf,['/shared/dondo/home/ml423/world_mts_runs/figs/synthesis/hima_AIRF_stats_' nts '.jpg'],'Resolution',700)
    end

    if cdfmap == 1
        locc = ismember(hima_AIstats(:,1),basinlist);
        figure
        t = tiledlayout(1,1,'Padding','none');
        t.Units = 'inches';
        t.OuterPosition = [0.25 0.25 4.95 3.6];
        nexttile;
        cd1 = cdfplot(hima_AIstats(locc==1,3));set(cd1,'Color','r','LineWidth',2);
        hold on
        cd2 = cdfplot(hima_AIstats(locc==1,4));set(cd2,'Color','b','LineWidth',2);
        hold off
        xlim([-3 1])
        ylim([0 1.02])
        xlabel('KGE')
        ylabel('CDF')
        title(' ')
        box on
        %legend(["ERA5L" "ERA5L-AI"],'Location','Northwest')
        set(gca,'FontName','Times New Roman','FontSize',14)
        exportgraphics(gcf,['/shared/dondo/home/ml423/world_mts_runs/figs/synthesis/KGE_CDF_comparison_AIRF_HIMA_' nts '.jpg'],'Resolution',600)
    end

end
