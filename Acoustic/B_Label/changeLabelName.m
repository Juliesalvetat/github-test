clear all
close all
clc;

% Load Label.mat 
load('D:\HAC_FAROFA_3\Cruise_FAROFA3\Treatment20200422_194552_FISH\Label\Label.mat')

path2Label = 'D:\HAC_FAROFA_3\Cruise_FAROFA3\Treatment20200624_150830_FISH\Label\Label.mat';
% path2Label = 'C:\Users\jsalveta\Desktop\Data\FAROFA_2\EK80\HAC\Cruise_FAROFA_2\Treatment20200120_010220_FISH\\Label\Label.mat';
% path2Label = 'C:\Users\jsalveta\Desktop\Data\FAROFA_1\EK80\HAC\Cruise_FAROFA_1\Treatment20191128_201456_FISH\Label\Label.mat';
load(path2Label);

LabelName = cell(14,1);


LabelName{1}='1-bottom.small.fish.school';
LabelName{2}='2-small.pelagics.school';
LabelName{3}='3-predators';
LabelName{4}='4-individual.demersal.fish';
LabelName{5}='5-shelfbreak.school';
LabelName{6}='6-shelfbreak.large.fish';
LabelName{7}='7-small.pelagics.and.predators';
LabelName{8}='8-mixt.reef.fish';
LabelName{9}='9-dense.school';
LabelName{10}='10-loose.school';
LabelName{11}='11-bottom.fish.thin.layer';
LabelName{12}='12-mel.nig';
LabelName{13}='13-can.suf';
LabelName{14}='14-Instrument';

save(path2Label,'LabelName','-append') 

save(path2Label,'LabelName','RemovedContour','ContourDepth','ContourIndPing','LabelIndex','-append') 

% 

% LabelName{1}='1-bottom.small.fish.school';
% LabelName{2}='2-sph.bar.vid';
% LabelName{3}='3-abu.sax.vid';
% LabelName{4}='4-mel.nig.vid';
% LabelName{5}='5-small.pelagics.school';
% LabelName{6}='6-Xstrong.echo.bottom';
% LabelName{7}='7-predators';
% LabelName{8}='8-shelfbreak.school';
% LabelName{9}='9-kyp.sec.vid';
% LabelName{10}='10-pelagic.predator';
% LabelName{11}='11-shelfbreak.large.fish';
% LabelName{12}='12-mel.nig';
% LabelName{13}='13-dec.mac.vid';
% LabelName{14}='14-Xmel_nig_or_caranx';
% LabelName{15}='15-ela.bip.vid';
% LabelName{16}='16-XUnknown.pelagic.fish';
% LabelName{17}='17-pelagic.predator.ela.bip.car.fal.vid';
% LabelName{18}='18-lut.joc.vid';
% LabelName{19}='19-car.rub.vid';
% LabelName{20}='20-individual.demersal.fish';
% LabelName{21}='21-mixt.reef.fish';
% LabelName{22}='22-small.pelagic.and.predators.vid';
% LabelName{23}='23-chr.mul.vid';
% LabelName{24}='24-damsel.grp';
% LabelName{25}='25-sand.fish';
% LabelName{26}='26-car.lat.vid';
% LabelName{27}='27-mel.nig.and.other.car.lat.sph.bar.car.cry.can.suf';
% LabelName{28}='28-mel.nig.lut.joc.sph.bar.vid';
% LabelName{29}='29-mel.nig.lut.joc.vid';
% LabelName{30}='30-can.suf.vid';
% LabelName{31}='31-mel.nig.car.cry.can.suf.vid';
% LabelName{32}='32-mel.nig.car.lat.sph.bar.aca.sol.vid';
% LabelName{33}='33-Instrument';
% LabelName{34}='34-Bootom.corr';
% LabelName{35}='35-bottom.small.fish';
% LabelName{36}='36-bottom.high200';
% LabelName{37}='37-layer.school';
% LabelName{38}='38-can.suf';
% LabelName{39}='39-car.lug.vid';
% LabelName{40}='40-can.suf.aca.sol.sph.bar.vid';
% LabelName{41}='41-car.cry.vid';
% LabelName{42}='42-can.suf.alu.scr.vid';
% LabelName{43}='43-ser.dum.gear';
% LabelName{44}='44-car.lug.gear';
% LabelName{45}='45-pelagic.predator.vid';
% LabelName{46}='46-bottom.fish.layer';
% LabelName{47}='47-aca.sol.vid';
% LabelName{48}='48-dense.school';
% LabelName{49}='49-loose.school';
% LabelName{50}='50-mel.nig.sph.bar.vid';
% LabelName{51}='51-small.pelagic.vid';
% LabelName{52}='52-bottom.school';
% LabelName{53}='53-small.pelagics.and.predators';
% LabelName{54}='54-unknow';
% LabelName{55}='55-y a du can.suf la dedans.vid';
% LabelName{56}='56-can.suf.mel.nig.alu.scr.car.cry.vid';
% LabelName{57}='57-bar.sph.gear';
% LabelName{58}='58-alu.scr.vid';
% LabelName{59}='59-car.lug.and.or.alu.scri';
% LabelName{60}='60-contains.can.suf.alu.scr.car.lug.ser.dum';
% LabelName{61}='61-can.suf.sph.bar.ela.bip.vid';

% Change indexes to regroupe labels
LabelIndex(LabelIndex==10)=11; %pelagic.predators become shelfbreak.large.fish
LabelIndex(LabelIndex==37|LabelIndex==25)=1;%blob
LabelIndex(LabelIndex==54|LabelIndex==35|LabelIndex==36)=46;%bottom.fish.layer
LabelIndex(LabelIndex==52)=21; %mixt.reef.fish
