clear all;
run=0;
%build animals script

%processing options (in this order as well)
params.bksub_each=1;
params.norm=1;
params.kalman=1;
params.plotrois=0;
params.imsize=[256 256];
params.nTrials=16;
params.framesPerTrial=60;
params.startFrame=9;
params.stimFrame=10;
params.trialTypes=11; %aka numel(MPa)

%order these according to the load order
MPa={0.2 0.4 0.6 0.8 0.8 0.8 0.8 0.8 0.8 0.8 0.8};
DC= {50  50  05  50  0   5   10  20  30  40  50};

mousedata=struct('name','','roidir','','maindir','','datadir','','MPa','','DC','','data','','params','');

%%
%%mouse 4
id=1;
mousedata(id).name='mouse4';
mousedata(id).roidir='/gpfs/data/shohamlab/shared_data/01012019/mouse4/RoiSet';
mousedata(id).maindatadir='/gpfs/data/shohamlab/shared_data/01012019/mouse4/';
mousedata(id).MPa=MPa;
mousedata(id).DC=DC;
mousedata(id).params=params;

%0.2MPa_50DC
mousedata(id).datadir{1}='TSeries-01132016-micea_2.1MHZ_3.18Hz_0.2MPa_50DC-250_mc';

%0.4MPa_50DC
mousedata(id).datadir{2}='TSeries-01132016-micea_2.1MHZ_3.18Hz_0.4MPa_50DC-249_mc';

%0.6MPa_50DC
mousedata(id).datadir{3}='TSeries-01132016-micea_2.1MHZ_3.18Hz_0.6MPa_50DC-248_mc';

%0.8MPa_0DC
mousedata(id).datadir{4}='TSeries-01132016-micea_2.1MHZ_3.18Hz_0.8MPa_0DC-247_mc';

%0.8MPa_5DC
mousedata(id).datadir{5}='TSeries-01132016-micea_2.1MHZ_3.18Hz_0.8MPa_5DC-246_mc';

%0.8MPa_10DC
mousedata(id).datadir{6}='TSeries-01132016-micea_2.1MHZ_3.18Hz_0.8MPa_10DC-245_mc';

%0.8MPa_20DC
mousedata(id).datadir{7}='TSeries-01132016-micea_2.1MHZ_3.18Hz_0.8MPa_20DC-244_mc';

%0.8MPa_30DC
mousedata(id).datadir{8}='TSeries-01132016-micea_2.1MHZ_3.18Hz_0.8MPa_30DC-243_mc';

%0.8MPa_40DC
mousedata(id).datadir{9}='TSeries-01132016-micea_2.1MHZ_3.18Hz_0.8MPa_40DC-242_mc';

%0.8MPa_50DC_1
mousedata(id).datadir{10}='TSeries-01132016-micea_2.1MHZ_3.18Hz_0.8MPa_50DC-241_mc';

%0.8MPa_50DC_2
%mousedata(id).datadir{11}='';


%%
%%mouse 5
id=2;
mousedata(id).name='mouse5';
mousedata(id).roidir='/gpfs/data/shohamlab/shared_data/01012019/mouse5/RoiSet';
mousedata(id).maindatadir='/gpfs/data/shohamlab/shared_data/01012019/mouse5';
mousedata(id).MPa=MPa;
mousedata(id).DC=DC;
mousedata(id).params=params;

%0.2MPa_50DC
mousedata(id).datadir{1}='TSeries-01132016-miceb_2.1MHZ_3.28Hz_0.2MPa_50DC-261_mc';

%0.4MPa_50DC
mousedata(id).datadir{2}='TSeries-01132016-miceb_2.1MHZ_3.28Hz_0.4MPa_50DC-260_mc';

%0.6MPa_50DC
mousedata(id).datadir{3}='TSeries-01132016-miceb_2.1MHZ_3.28Hz_0.6MPa_50DC-259_mc';

%0.8MPa_0DC
mousedata(id).datadir{4}='TSeries-01132016-miceb_2.1MHZ_3.28Hz_0.8MPa_0DC-258_mc';

%0.8MPa_5DC
mousedata(id).datadir{5}='TSeries-01132016-miceb_2.1MHZ_3.28Hz_0.8MPa_5DC-257_mc';

%0.8MPa_10DC
mousedata(id).datadir{6}='TSeries-01132016-miceb_2.1MHZ_3.28Hz_0.8MPa_10DC-256_mc';

%0.8MPa_20DC
mousedata(id).datadir{7}='TSeries-01132016-miceb_2.1MHZ_3.28Hz_0.8MPa_20DC-255_mc';

%0.8MPa_30DC
mousedata(id).datadir{8}='TSeries-01132016-miceb_2.1MHZ_3.28Hz_0.8MPa_30DC-254_mc';

%0.8MPa_40DC
mousedata(id).datadir{9}='TSeries-01132016-miceb_2.1MHZ_3.28Hz_0.8MPa_40DC-253_mc';

%0.8MPa_50DC_1
mousedata(id).datadir{10}='TSeries-01132016-miceb_2.1MHZ_3.28Hz_0.8MPa_50DC-252_mc';

%0.8MPa_50DC_2
%mousedata(id).datadir{11}='';

%%
%%mouse 7
id=3;
mousedata(id).name='mouse7';
mousedata(id).roidir='/gpfs/data/shohamlab/shared_data/01012019/mouse7/RoiSet';
mousedata(id).maindatadir='/gpfs/data/shohamlab/shared_data/01012019/mouse7';
mousedata(id).MPa=MPa;
mousedata(id).DC=DC;
mousedata(id).params=params;

%0.2MPa_50DC
mousedata(id).datadir{1}='TSeries-01132016-miceb_2.1MHz_3.18Hz_0.2MPa_50DC-294_mc';

%0.4MPa_50DC
mousedata(id).datadir{2}='TSeries-01132016-miceb_2.1MHz_3.18Hz_0.4MPa_50DC-293_mc';

%0.6MPa_50DC
mousedata(id).datadir{3}='TSeries-01132016-miceb_2.1MHz_3.18Hz_0.6MPa_50DC-292_mc';

%0.8MPa_0DC
mousedata(id).datadir{4}='TSeries-01132016-miceb_2.1MHz_3.18Hz_0.8MPa_0DC-290_mc';

%0.8MPa_5DC
mousedata(id).datadir{5}='TSeries-01132016-miceb_2.1MHz_3.18Hz_0.8MPa_5DC-289_mc';

%0.8MPa_10DC
mousedata(id).datadir{6}='TSeries-01132016-miceb_2.1MHz_3.18Hz_0.8MPa_10DC-288_mc';

%0.8MPa_20DC
mousedata(id).datadir{7}='TSeries-01132016-miceb_2.1MHz_3.18Hz_0.8MPa_20DC-287_mc';

%0.8MPa_30DC
mousedata(id).datadir{8}='TSeries-01132016-miceb_2.1MHz_3.18Hz_0.8MPa_30DC-286_mc';

%0.8MPa_40DC
mousedata(id).datadir{9}='TSeries-01132016-miceb_2.1MHz_3.18Hz_0.8MPa_40DC-285_mc';

%0.8MPa_50DC_1
%datadir{10}='TSeries-01132016-miceb_2.1MHz_3.18Hz_0.8MPa_50DC-283_mc';
mousedata(id).datadir{10}='TSeries-01132016-miceb_2.1MHz_3.18Hz_0.8MPa_50DC-284_mc';
%0.8MPa_50DC_2
mousedata(id).datadir{11}='TSeries-01132016-miceb_2.1MHz_3.18Hz_0.8MPa_50DC-291_mc';

%%
%mouse 8
id=4;
mousedata(id).name='mouse8';
mousedata(id).roidir='/gpfs/data/shohamlab/shared_data/01012019/mouse8/RoiSet';
mousedata(id).maindatadir='/gpfs/data/shohamlab/shared_data/01012019/mouse8';
mousedata(id).MPa=MPa;
mousedata(id).DC=DC;
mousedata(id).params=params;

%0.2MPa_50DC
mousedata(id).datadir{1}='TSeries-01132016-miced_3.46Hz_2.1MHz_0.2MPa_50DC-287_mc';

%0.4MPa_50DC
mousedata(id).datadir{2}='TSeries-01132016-miced_3.46Hz_2.1MHz_0.4MPa_50DC-286_mc';

%0.6MPa_50DC
mousedata(id).datadir{3}='TSeries-01132016-miced_3.46Hz_2.1MHz_0.6MPa_50DC-285_mc';

%0.8MPa_0DC
mousedata(id).datadir{4}='TSeries-01132016-miced_3.46Hz_2.1MHz_0.8MPa_0DC-283_mc';

%0.8MPa_5DC
mousedata(id).datadir{5}='TSeries-01132016-miced_3.46Hz_2.1MHz_0.8MPa_5DC-282_mc';

%0.8MPa_10DC
mousedata(id).datadir{6}='TSeries-01132016-miced_3.46Hz_2.1MHz_0.8MPa_10DC-281_mc';

%0.8MPa_20DC
mousedata(id).datadir{7}='TSeries-01132016-miced_3.46Hz_2.1MHz_0.8MPa_20DC-280_mc';

%0.8MPa_30DC
mousedata(id).datadir{8}='TSeries-01132016-miced_3.46Hz_2.1MHz_0.8MPa_30DC-279_mc';

%0.8MPa_40DC
mousedata(id).datadir{9}='TSeries-01132016-miced_3.46Hz_2.1MHz_0.8MPa_40DC-278_mc';

%0.8MPa_50DC_1
mousedata(id).datadir{10}='TSeries-01132016-miced_3.46Hz_2.1MHz_0.8MPa_50DC-277_mc';

%0.8MPa_50DC_2
%mousedata(id).datadir{11}='';

%%
%mouse 9
id=5;
mousedata(id).name='mouse9';
mousedata(id).roidir='/gpfs/data/shohamlab/shared_data/01012019/mouse9/RoiSet';
mousedata(id).maindatadir='/gpfs/data/shohamlab/shared_data/01012019/mouse9';
mousedata(id).MPa=MPa;
mousedata(id).DC=DC;
mousedata(id).params=params;

%0.2MPa_50DC
mousedata(id).datadir{1}='TSeries-01132016-micec_3.46Hz_2.1MHz_0.2MPa_50DC-275_mc';

%0.4MPa_50DC
mousedata(id).datadir{2}='TSeries-01132016-micec_3.46Hz_2.1MHz_0.4MPa_50DC-274_mc';

%0.6MPa_50DC
mousedata(id).datadir{3}='TSeries-01132016-micec_3.46Hz_2.1MHz_0.6MPa_50DC-273_mc';

%0.8MPa_0DC
mousedata(id).datadir{4}='TSeries-01132016-micec_3.46Hz_2.1MHz_0.8MPa_0DC-271_mc';

%0.8MPa_5DC
mousedata(id).datadir{5}='TSeries-01132016-micec_3.46Hz_2.1MHz_0.8MPa_5DC-270_mc';

%0.8MPa_10DC
mousedata(id).datadir{6}='TSeries-01132016-micec_3.46Hz_2.1MHz_0.8MPa_10DC-269_mc';

%0.8MPa_20DC
mousedata(id).datadir{7}='TSeries-01132016-micec_3.46Hz_2.1MHz_0.8MPa_20DC-268_mc';

%0.8MPa_30DC
mousedata(id).datadir{8}='TSeries-01132016-micec_3.46Hz_2.1MHz_0.8MPa_30DC-267_mc';

%0.8MPa_40DC
mousedata(id).datadir{9}='TSeries-01132016-micec_3.46Hz_2.1MHz_0.8MPa_40DC-266_mc';

%0.8MPa_50DC_1
mousedata(id).datadir{10}='TSeries-01132016-micec_3.46Hz_2.1MHz_0.8MPa_50DC-272_mc';

%0.8MPa_50DC_2
%datadir{11}='TSeries-01132016-micec_2.1MHz_3.375HZ_0.8MPa_50DC-265_mc';
%mousedata(id).datadir{11}='';

%%
%%mouse 10
id=6;
mousedata(id).name='mouse10';
mousedata(id).roidir='/gpfs/data/shohamlab/shared_data/01012019/mouse10/RoiSet';
mousedata(id).maindatadir='/gpfs/data/shohamlab/shared_data/01012019/mouse10';
mousedata(id).MPa=MPa;
mousedata(id).DC=DC;
mousedata(id).params=params;

%0.2MPa_50DC
mousedata(id).datadir{1}='TSeries-01132016-miceb_2.1MHz_3.375Hz_0.2MPa_50DC-261_mc';

%0.4MPa_50DC
mousedata(id).datadir{2}='TSeries-01132016-miceb_2.1MHz_3.375Hz_0.4MPa_50DC-260_mc';

%0.6MPa_50DC
mousedata(id).datadir{3}='TSeries-01132016-miceb_2.1MHz_3.375Hz_0.6MPa_50DC-259_mc';

%0.8MPa_0DC
mousedata(id).datadir{4}='TSeries-01132016-miceb_2.1MHz_3.375HZ_0.8MPa_0DC-257_mc';

%0.8MPa_5DC
mousedata(id).datadir{5}='TSeries-01132016-miceb_2.1MHz_3.375HZ_0.8MPa_5DC-256_mc';

%0.8MPa_10DC
mousedata(id).datadir{6}='TSeries-01132016-miceb_2.1MHz_3.375HZ_0.8MPa_10DC-255_mc';

%0.8MPa_20DC
mousedata(id).datadir{7}='TSeries-01132016-miceb_2.1MHz_3.375HZ_0.8MPa_20DC-254_mc';

%0.8MPa_30DC
mousedata(id).datadir{8}='TSeries-01132016-miceb_2.1MHz_3.375HZ_0.8MPa_30DC-253_mc';

%0.8MPa_40DC
mousedata(id).datadir{9}='TSeries-01132016-miceb_2.1MHz_3.375HZ_0.8MPa_40DC-252_mc';

%0.8MPa_50DC_1
mousedata(id).datadir{10}='TSeries-01132016-miceb_2.1MHz_3.375HZ_0.8MPa_50DC-251_mc';

%0.8MPa_50DC_2
mousedata(id).datadir{11}='TSeries-01132016-miceb_2.1MHz_3.375HZ_0.8MPa_50DC-258_mc';

%%
%%mouse 11
id=7;
mousedata(id).name='mouse11';
mousedata(id).roidir='/gpfs/data/shohamlab/shared_data/01012019/mouse11/RoiSet';
mousedata(id).maindatadir='/gpfs/data/shohamlab/shared_data/01012019/mouse11';
mousedata(id).MPa=MPa;
mousedata(id).DC=DC;
mousedata(id).params=params;

%0.2MPa_50DC
mousedata(id).datadir{1}='TSeries-01132016-micea_2.1MHz_3.375HZ_0.2MPa_50DC-248_mc';

%0.4MPa_50DC
mousedata(id).datadir{2}='TSeries-01132016-micea_2.1MHz_3.375HZ_0.4MPa_50DC-247_mc';

%0.6MPa_50DC
mousedata(id).datadir{3}='TSeries-01132016-micea_2.1MHz_3.375HZ_0.6MPa_50DC-246_mc';

%0.8MPa_0DC
mousedata(id).datadir{4}='TSeries-01132016-micea_2.1MHz_3.375HZ_0.8MPa_0DC-244_mc';

%0.8MPa_5DC
mousedata(id).datadir{5}='TSeries-01132016-micea_2.1MHz_3.375HZ_0.8MPa_5DC-243_mc';

%0.8MPa_10DC
mousedata(id).datadir{6}='TSeries-01132016-micea_2.1MHz_3.375HZ_0.8MPa_10DC-242_mc';

%0.8MPa_20DC
mousedata(id).datadir{7}='TSeries-01132016-micea_2.1MHz_3.375HZ_0.8MPa_20DC-241_mc';

%0.8MPa_30DC
mousedata(id).datadir{8}='TSeries-01132016-micea_2.1MHz_3.375HZ_0.8MPa_30DC-240_mc';

%0.8MPa_40DC
mousedata(id).datadir{9}='TSeries-01132016-micea_2.1MHz_3.375HZ_0.8MPa_40DC-239_mc';

%0.8MPa_50DC_1
mousedata(id).datadir{10}='TSeries-01132016-micea_2.1MHz_3.375HZ_0.8MPa_50DC-238_mc';

%0.8MPa_50DC_1
mousedata(id).datadir{11}='TSeries-01132016-micea_2.1MHz_3.375HZ_0.8MPa_50DC-245_mc';
%%



%Run analysis
if run==1
yi_analysis(mousedata)
end



