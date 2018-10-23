%% TOFマップ作成 Group1
clear; close all;
for i = 1:5
    cd('H:\result\2018_02_28_-kwave');
    myfilename = sprintf('2018_04_26_Group1ERate5IRate02Num%d', i);
    cd(myfilename);
    load('kgrid.mat')
    load('rfdata.mat')
    tof_data = get_tof_AIC_from_singleRFData(rfdata,kgrid,75);
    cd('H:\result\2018_04_27_analyzeLipidModel')
    myfilename = sprintf('2018_04_27_TOFdata_Group1ERate5IRate02Num%d', i);
    save(myfilename,'tof_data');
    clear; close all;
end
for i = 1:4
    cd('H:\result\2018_02_28_-kwave');
    myfilename = sprintf('2018_04_26_Group1ERate%dIRate02Num3', i);
    cd(myfilename);
    load('kgrid.mat')
    load('rfdata.mat')
    tof_data = get_tof_AIC_from_singleRFData(rfdata,kgrid,75);
    cd('H:\result\2018_04_27_analyzeLipidModel')
    myfilename = sprintf('2018_04_27_TOFdata_Group1ERate%dIRate02Num3', i);
    save(myfilename,'tof_data');
    clear; close all;
end
for i = 4:2:8
    cd('H:\result\2018_02_28_-kwave');
    myfilename = sprintf('2018_04_26_Group1ERate5IRate0%dNum3', i);
    cd(myfilename);
    load('kgrid.mat')
    load('rfdata.mat')
    tof_data = get_tof_AIC_from_singleRFData(rfdata,kgrid,75);
    cd('H:\result\2018_04_27_analyzeLipidModel')
    myfilename = sprintf('2018_04_27_TOFdata_Group1ERate5IRate0%dNum3', i);
    save(myfilename,'tof_data');
    clear; close all;
end

cd('H:\result\2018_02_28_-kwave');
myfilename = sprintf('2018_04_26_Group1ERate5IRate10Num3');
cd(myfilename);
load('kgrid.mat')
load('rfdata.mat')
tof_data = get_tof_AIC_from_singleRFData(rfdata,kgrid,75);
cd('H:\result\2018_04_27_analyzeLipidModel')
myfilename = sprintf('2018_04_27_TOFdata_Group1ERate5IRate10Num3');
save(myfilename,'tof_data');
clear; close all;

%% 平均音速算出(配列定義) Group1
clear; close all;
tof_cell  = zeros(200,100,5,5,5);%{Receiver, Transmitter, ERate, IRate, Num, }
aveSOS = zeros(5,5,5);%{ERate, IRate, Num}
for i = 1:5
    cd('H:\result\2018_04_27_analyzeLipidModel')
    myfilename = sprintf('2018_04_27_TOFdata_Group1ERate5IRate02Num%d', i);
    load(myfilename,'tof_data');
    tof_cell(:,:,5,1,i) = tof_data;
end
for i = 1:4
    cd('H:\result\2018_04_27_analyzeLipidModel')
    myfilename = sprintf('2018_04_27_TOFdata_Group1ERate%dIRate02Num3', i);
    load(myfilename,'tof_data');
    tof_cell(:,:,i,1,3) = tof_data;
end
for i = 4:2:8
    cd('H:\result\2018_04_27_analyzeLipidModel')
    myfilename = sprintf('2018_04_27_TOFdata_Group1ERate5IRate0%dNum3', i);
    load(myfilename,'tof_data');
    tof_cell(:,:,5,i/2,3) = tof_data;
end
cd('H:\result\2018_04_27_analyzeLipidModel')
myfilename = sprintf('2018_04_27_TOFdata_Group1ERate5IRate10Num3');
load(myfilename,'tof_data');
tof_cell(:,:,5,5,3) = tof_data;
myfilename = sprintf('2018_04_27_Group1allTOFdata_matrix');
save(myfilename,'tof_cell');
%% 平均音速算出（メイン） Group1
cd('H:\result\2018_04_27_analyzeLipidModel')
myfilename = sprintf('2018_04_27_Group1allTOFdata_matrix');
load(myfilename,'tof_cell');
cd('H:\result\2018_02_28_-kwave');
myfilename = sprintf('2018_04_26_Group1ERate5IRate02Num3');
cd(myfilename);
load('kgrid.mat')
param = makeParam( 0.125,0.125,400,400,0 );
param.source.point_map = cast(linspace(1,100,100),'int8');
t_size = 40.e-3;
t_num = 200;%トランスデューサ数
t_pos = zeros(2, t_num);%センサ位置
t_pos(1,1:t_num/2) = -t_size/2:t_size/(t_num/2-1):t_size/2 ;%素子水平方向距離[m]
t_pos(2,1:t_num/2) = t_size/2;
t_pos(1,t_num/2+1:t_num) = t_pos(1,1:t_num/2);
t_pos(2,t_num/2+1:t_num) = -t_size/2;
sensor.mask=t_pos;
leng = zeros(1,100);
for i = 1:5
    [Min,ind]= min(tof_cell(101:end,:,i,1,3));
    ind = ind+100;
    for j = 1: t_num/2
        leng(1,j) = norm(t_pos(:,j)-t_pos(:,ind(1,j)));
    end
    aveSOS(i,1,3) = mean(leng./Min);
end
for i = 1:5
    if i == 3
        continue
    end
    [Min,ind]= min(tof_cell(101:end,:,5,1,i));
    ind = ind+100;
    for j = 1: t_num/2
        leng(1,j) = norm(t_pos(:,j)-t_pos(:,ind(1,j)));
    end
    aveSOS(5,1,i) = mean(leng./Min);
end
for i = 2:5
    [Min,ind]= min(tof_cell(101:end,:,5,i,3));
    ind = ind+100;
    for j = 1: t_num/2
        leng(1,j) = norm(t_pos(:,j)-t_pos(:,ind(1,j)));
    end
    aveSOS(5,i,3) = mean(leng./Min);
end
cd('H:\result\2018_04_27_analyzeLipidModel')
myfilename = sprintf('2018_04_28_Group1aveSOS');
save(myfilename,'aveSOS');

%% TOFマップ作成
%Group2
clear; close all;
for i = 1:5
    cd('H:\result\2018_02_28_-kwave');
    myfilename = sprintf('2018_04_26_Group2ERate5IRate02Num%d', i);
    cd(myfilename);
    load('kgrid.mat')
    load('rfdata.mat')
    tof_data = get_tof_AIC_from_singleRFData(rfdata,kgrid,75);
    cd('H:\result\2018_04_27_analyzeLipidModel')
    myfilename = sprintf('2018_04_27_TOFdata_Group2ERate5IRate02Num%d', i);
    save(myfilename,'tof_data');
    clear; close all;
end
for i = 1:4
    cd('H:\result\2018_02_28_-kwave');
    myfilename = sprintf('2018_04_26_Group2ERate%dIRate02Num3', i);
    cd(myfilename);
    load('kgrid.mat')
    load('rfdata.mat')
    tof_data = get_tof_AIC_from_singleRFData(rfdata,kgrid,75);
    cd('H:\result\2018_04_27_analyzeLipidModel')
    myfilename = sprintf('2018_04_27_TOFdata_Group2ERate%dIRate02Num3', i);
    save(myfilename,'tof_data');
    clear; close all;
end
for i = 4:2:8
    cd('H:\result\2018_02_28_-kwave');
    myfilename = sprintf('2018_04_26_Group2ERate5IRate0%dNum3', i);
    cd(myfilename);
    load('kgrid.mat')
    load('rfdata.mat')
    tof_data = get_tof_AIC_from_singleRFData(rfdata,kgrid,75);
    cd('H:\result\2018_04_27_analyzeLipidModel')
    myfilename = sprintf('2018_04_27_TOFdata_Group2ERate5IRate0%dNum3', i);
    save(myfilename,'tof_data');
    clear; close all;
end

cd('H:\result\2018_02_28_-kwave');
myfilename = sprintf('2018_04_26_Group2ERate5IRate10Num3');
cd(myfilename);
load('kgrid.mat')
load('rfdata.mat')
tof_data = get_tof_AIC_from_singleRFData(rfdata,kgrid,75);
cd('H:\result\2018_04_27_analyzeLipidModel')
myfilename = sprintf('2018_04_27_TOFdata_Group2ERate5IRate10Num3');
save(myfilename,'tof_data');
clear; close all;

%% 平均音速算出(配列定義) Group2
clear; close all;
tof_cell  = zeros(200,100,5,5,5);%{Receiver, Transmitter, ERate, IRate, Num, }
aveSOS = zeros(5,5,5);%{ERate, IRate, Num}
for i = 1:5
    cd('H:\result\2018_04_27_analyzeLipidModel')
    myfilename = sprintf('2018_04_27_TOFdata_Group2ERate5IRate02Num%d', i);
    load(myfilename,'tof_data');
    tof_cell(:,:,5,1,i) = tof_data;
end
for i = 1:4
    cd('H:\result\2018_04_27_analyzeLipidModel')
    myfilename = sprintf('2018_04_27_TOFdata_Group2ERate%dIRate02Num3', i);
    load(myfilename,'tof_data');
    tof_cell(:,:,i,1,3) = tof_data;
end
for i = 4:2:8
    cd('H:\result\2018_04_27_analyzeLipidModel')
    myfilename = sprintf('2018_04_27_TOFdata_Group2ERate5IRate0%dNum3', i);
    load(myfilename,'tof_data');
    tof_cell(:,:,5,i/2,3) = tof_data;
end
cd('H:\result\2018_04_27_analyzeLipidModel')
myfilename = sprintf('2018_04_27_TOFdata_Group2ERate5IRate10Num3');
load(myfilename,'tof_data');
tof_cell(:,:,5,5,3) = tof_data;
myfilename = sprintf('2018_04_27_Group2allTOFdata_matrix');
save(myfilename,'tof_cell');

%% 平均音速算出（メイン）Group2
cd('H:\result\2018_04_27_analyzeLipidModel')
myfilename = sprintf('2018_04_27_Group2allTOFdata_matrix');
load(myfilename,'tof_cell');
cd('H:\result\2018_02_28_-kwave');
myfilename = sprintf('2018_04_26_Group2ERate5IRate02Num3');
cd(myfilename);
load('kgrid.mat')
param = makeParam( 0.125,0.125,400,400,0 );
param.source.point_map = cast(linspace(1,100,100),'int8');
t_size = 40.e-3;
t_num = 200;%トランスデューサ数
t_pos = zeros(2, t_num);%センサ位置
t_pos(1,1:t_num/2) = -t_size/2:t_size/(t_num/2-1):t_size/2 ;%素子水平方向距離[m]
t_pos(2,1:t_num/2) = t_size/2;
t_pos(1,t_num/2+1:t_num) = t_pos(1,1:t_num/2);
t_pos(2,t_num/2+1:t_num) = -t_size/2;
sensor.mask=t_pos;
leng = zeros(1,100);
for i = 1:5
    [Min,ind]= min(tof_cell(101:end,:,i,1,3));
    ind = ind+100;
    for j = 1: t_num/2
        leng(1,j) = norm(t_pos(:,j)-t_pos(:,ind(1,j)));
    end
    aveSOS(i,1,3) = mean(leng./Min);
end
for i = 1:5
    if i == 3
        continue
    end
    [Min,ind]= min(tof_cell(101:end,:,5,1,i));
    ind = ind+100;
    for j = 1: t_num/2
        leng(1,j) = norm(t_pos(:,j)-t_pos(:,ind(1,j)));
    end
    aveSOS(5,1,i) = mean(leng./Min);
end
for i = 2:5
    [Min,ind]= min(tof_cell(101:end,:,5,i,3));
    ind = ind+100;
    for j = 1: t_num/2
        leng(1,j) = norm(t_pos(:,j)-t_pos(:,ind(1,j)));
    end
    aveSOS(5,i,3) = mean(leng./Min);
end
cd('H:\result\2018_04_27_analyzeLipidModel')
myfilename = sprintf('2018_04_28_Group2aveSOS');
save(myfilename,'aveSOS');

%% Group1,Group2のaveSOSを平均してみる．
%ねらい：結果の非線形さがモデルのランダムさによらず同じ傾向を示すのならば，
%             算出アルゴリズムになにか誤りがあるのではないかと疑えるので検証する
%計算
cd('H:\result\2018_04_27_analyzeLipidModel')
load('2018_04_28_Group1aveSOS.mat')
aveSOS_Gr1 = aveSOS;
load('2018_04_28_Group2aveSOS.mat')
aveSOS_Gr2 = aveSOS;
aveSOS_Gr1_Gr2 = (aveSOS_Gr1 + aveSOS_Gr2) /2;
%グラフ作成
figure;
plot(aveSOS_Gr1_Gr2(:,1,3));
%% make figure置き場　[2018-04-30]
aveSOS_num =  zeros(1,5);
aveSOS_num(1,:) = aveSOS(5,1,:);
figure;
plot(aveSOS_num);
xlabel('EMCLの個数');
ylabel('平均音速[m/s]')
exportfig('H:\result\2018_04_27_analyzeLipidModel\Group1aveSOS_num','png',[300,200]);

aveSOS_EMCLrate = zeros(1,5);
aveSOS_EMCLrate(1,:) = aveSOS(:,1,3);
figure;
plot(aveSOS_EMCLrate);
xlabel('EMCLの占有率[%]');
ylabel('平均音速[m/s]')
exportfig('H:\result\2018_04_27_analyzeLipidModel\Group1aveSOS_EMCLrate','png',[300,200]);

aveSOS_IMCLrate = zeros(1,5);
aveSOS_IMCLrate(1,:) = aveSOS(5,:,3);
figure;
x = [ 0.1 0.2 0.3 0.4 0.5 ];
plot(x,aveSOS_IMCLrate);
xlabel('IMCLの占有率[%]');
ylabel('平均音速[m/s]')
exportfig('H:\result\2018_04_27_analyzeLipidModel\Group1aveSOS_IMCLrate','png',[300,200]);

aveSOS_num =  zeros(1,5);
aveSOS_num(1,:) = aveSOS(5,1,:);
figure;
plot(aveSOS_num);
xlabel('EMCLの個数');
ylabel('平均音速[m/s]')
exportfig('H:\result\2018_04_27_analyzeLipidModel\Group2aveSOS_num','png',[300,200]);
aveSOS_EMCLrate = zeros(1,5);
aveSOS_EMCLrate(1,:) = aveSOS(:,1,3);
figure;
plot(aveSOS_EMCLrate);
xlabel('EMCLの占有率[%]');
ylabel('平均音速[m/s]')
exportfig('H:\result\2018_04_27_analyzeLipidModel\Group2aveSOS_EMCLrate','png',[300,200]);

aveSOS_IMCLrate = zeros(1,5);
aveSOS_IMCLrate(1,:) = aveSOS(5,:,3);
figure;
x = [ 0.1 0.2 0.3 0.4 0.5 ];
plot(x,aveSOS_IMCLrate);
xlabel('IMCLの占有率[%]');
ylabel('平均音速[m/s]')
exportfig('H:\result\2018_04_27_analyzeLipidModel\Group2aveSOS_IMCLrate','png',[300,200]);


aveSOS_num =  zeros(1,5);
aveSOS_num(1,:) = aveSOS(5,1,:);
figure;
plot(aveSOS_num);
xlabel('EMCLの個数');
ylabel('平均音速[m/s]')
exportfig('H:\result\2018_04_27_analyzeLipidModel\Group1&2aveSOS_num','png',[300,200]);
aveSOS_EMCLrate = zeros(1,5);
aveSOS_EMCLrate(1,:) = aveSOS(:,1,3);
figure;
plot(aveSOS_EMCLrate);
xlabel('EMCLの占有率[%]');
ylabel('平均音速[m/s]')
exportfig('H:\result\2018_04_27_analyzeLipidModel\Group1&2aveSOS_EMCLrate','png',[300,200]);

aveSOS_IMCLrate = zeros(1,5);
aveSOS_IMCLrate(1,:) = aveSOS(5,:,3);
figure;
x = [ 0.1 0.2 0.3 0.4 0.5 ];
plot(x,aveSOS_IMCLrate);
xlabel('IMCLの占有率[%]');
ylabel('平均音速[m/s]')
exportfig('H:\result\2018_04_27_analyzeLipidModel\Group1&2aveSOS_IMCLrate','png',[300,200]);

%% EMCLの占有率を1,2,3,4,5,6,8,10,12,14,16,18,20%と振った場合
aveSOS_EMCLrate = zeros(1,20);
aveSOS_EMCLrate(1,:) = aveSOS(:,1,3);
figure;
x = [2 4 6 8 10 12 14 16 18 20];
plot(x(1,1:end),aveSOS_EMCLrate(1,2:2:20));
xlabel('EMCLの占有率[%]');
ylabel('平均音速[m/s]')
exportfig('H:\result\2018_04_27_analyzeLipidModel\Group1aveSOS_expand_EMCLrate','png',[300,200]);

%% 可視化段階
%subplotを想定．mediumとTOFマップ（最短到達時間ポイントに大きな値を埋め込んで図を見やすくする役割）
% 対象はGroup1に限定

param = makeParam( 0.125,0.125,400,400,0 );
kgrid = makeGrid(param.grid.Nx, param.grid.dx, param.grid.Ny, param.grid.dy);
for i = 1:5
    figure;
    cd('H:\data\kwave\medium\2018_04_26_randomScatter')
    myfilename = sprintf('Group1ERate%dIRate02Num3', i);
    load(myfilename);
    subplot(1,2,1);
    imagesc(kgrid.x_vec*1000,kgrid.y_vec*1000,medium.sound_speed);
    axis equal
    axis tight
    hold on
    scatter(t_pos(2,:)*1000,t_pos(1,:)*1000,'rs','filled');
    xlabel('y-axis[mm]')
    ylabel('x-axis[mm]')
    hold off
    subplot(1,2,2);
    cd('H:\result\2018_04_27_analyzeLipidModel')
    myfilename = sprintf('2018_04_27_TOFdata_Group1ERate%dIRate02Num3',i);
    load(myfilename)
    [Min,ind]= min(tof_data(101:end,:));
    ind = ind+100;
    for j = 1:100
        tof_data(ind(1,j),j) = 0;
    end
    imagesc(1:100,1:100,tof_data(101:end,:));
    xlabel('receiver')
    ylabel('transimitter')
    axis tight
    axis equal
    titlename = sprintf('Group1ERate%dIRate02Num3',i);
    axes;
    title(titlename);
    axis off;
    cd('H:\result\2018_04_27_analyzeLipidModel')
    myfilename = sprintf('2018_05_01_Group1ERate%dIRate02Num3_Compare',i);
    exportfig(myfilename,'png',[400,180]);
end

for i = 1:5
    figure;
    cd('H:\data\kwave\medium\2018_04_26_randomScatter')
    myfilename = sprintf('Group1ERate5IRate02Num%d', i);
    load(myfilename);
    subplot(1,2,1);
    imagesc(kgrid.x_vec*1000,kgrid.y_vec*1000,medium.sound_speed);
    axis equal
    axis tight
    hold on
    scatter(t_pos(2,:)*1000,t_pos(1,:)*1000,'rs','filled');
    xlabel('y-axis[mm]')
    ylabel('x-axis[mm]')
    hold off
    subplot(1,2,2);
    cd('H:\result\2018_04_27_analyzeLipidModel')
    myfilename = sprintf('2018_04_27_TOFdata_Group1ERate5IRate02Num%d',i);
    load(myfilename)
    [Min,ind]= min(tof_data(101:end,:));
    ind = ind+100;
    for j = 1:100
        tof_data(ind(1,j),j) = 0;
    end
    imagesc(1:100,1:100,tof_data(101:end,:));
    xlabel('receiver')
    ylabel('transimitter')
    axis tight
    axis equal
    titlename = sprintf('Group1ERate5IRate02Num%d',i);
    axes;
    title(titlename);
    axis off;
    cd('H:\result\2018_04_27_analyzeLipidModel')
    myfilename = sprintf('2018_05_01_Group1ERate5IRate02Num%d_Compare',i);
    exportfig(myfilename,'png',[400,180]);
end

for i = 1:5
    figure;
    cd('H:\data\kwave\medium\2018_04_26_randomScatter')
    myfilename = sprintf('Group1ERate5IRate%02dNum3', 2*i);
    load(myfilename);
    subplot(1,2,1);
    imagesc(kgrid.x_vec*1000,kgrid.y_vec*1000,medium.sound_speed);
    axis equal
    axis tight
    hold on
    scatter(t_pos(2,:)*1000,t_pos(1,:)*1000,'rs','filled');
    xlabel('y-axis[mm]')
    ylabel('x-axis[mm]')
    hold off
    subplot(1,2,2);
    cd('H:\result\2018_04_27_analyzeLipidModel')
    myfilename = sprintf('2018_04_27_TOFdata_Group1ERate5IRate%02dNum3',i*2);
    load(myfilename)
    [Min,ind]= min(tof_data(101:end,:));
    ind = ind+100;
    for j = 1:100
        tof_data(ind(1,j),j) = 0;
    end
    imagesc(1:100,1:100,tof_data(101:end,:));
    xlabel('receiver')
    ylabel('transimitter')
    axis tight
    axis equal
    titlename = sprintf('Group1ERate5IRate%02dNum3',i*2);
    axes;
    title(titlename);
    axis off;
    cd('H:\result\2018_04_27_analyzeLipidModel')
    myfilename = sprintf('2018_05_01_Group1ERate5IRate%02dNum3_Compare',i);
    exportfig(myfilename,'png',[400,180]);
end

for i = 2:4:20
    figure;
    cd('H:\data\kwave\medium\2018_04_26_randomScatter')
    myfilename = sprintf('Group1ERate%dIRate02Num3', i);
    load(myfilename);
    subplot(1,2,1);
    imagesc(kgrid.x_vec*1000,kgrid.y_vec*1000,medium.sound_speed);
    axis equal
    axis tight
    hold on
    scatter(t_pos(2,:)*1000,t_pos(1,:)*1000,'rs','filled');
    xlabel('y-axis[mm]')
    ylabel('x-axis[mm]')
    hold off
    subplot(1,2,2);
    cd('H:\result\2018_04_27_analyzeLipidModel')
    myfilename = sprintf('2018_04_27_TOFdata_Group1ERate%dIRate02Num3',i);
    load(myfilename)
    [Min,ind]= min(tof_data(101:end,:));
    ind = ind+100;
    for j = 1:100
        tof_data(ind(1,j),j) = 0;
    end
    imagesc(1:100,1:100,tof_data(101:end,:));
    xlabel('receiver')
    ylabel('transimitter')
    axis tight
    axis equal
    titlename = sprintf('Group1ERate%dIRate02Num3',i);
    axes;
    title(titlename);
    axis off;
    cd('H:\result\2018_04_27_analyzeLipidModel')
    myfilename = sprintf('2018_05_01_Group1ERate%dIRate02Num3_Compare',i);
    exportfig(myfilename,'png',[400,180]);
end

%% Mediumの作り直し
medium = makeRandomMedium(param,20,3,2);
imagesc(kgrid.x_vec*1000,kgrid.y_vec*1000,medium.sound_speed);
axis equal
axis tight
hold on
scatter(t_pos(2,:)*1000,t_pos(1,:)*1000,'rs','filled');
xlabel('y-axis[mm]')
ylabel('x-axis[mm]')
hold off
%% 可視化段階()
%subplotを想定．mediumとTOFマップ（最短到達時間ポイントに大きな値を埋め込んで図を見やすくする役割）
% 対象はGroup1に限定

param = makeParam( 0.125,0.125,400,400,0 );
kgrid = makeGrid(param.grid.Nx, param.grid.dx, param.grid.Ny, param.grid.dy);
for i = 1:5
    figure;
    cd('H:\data\kwave\medium\2018_04_26_randomScatter')
    myfilename = sprintf('Group1ERate%dIRate02Num3-2', i);
    load(myfilename);
    subplot(1,2,1);
    imagesc(kgrid.x_vec*1000,kgrid.y_vec*1000,medium.sound_speed);
    axis equal
    axis tight
    hold on
    scatter(t_pos(2,:)*1000,t_pos(1,:)*1000,'rs','filled');
    xlabel('y-axis[mm]')
    ylabel('x-axis[mm]')
    hold off
    subplot(1,2,2);
    cd('H:\result\2018_04_27_analyzeLipidModel')
    myfilename = sprintf('2018_04_27_TOFdata_Group1ERate%dIRate02Num3-2',i);
    load(myfilename)
    [Min,ind]= min(tof_data(101:end,:));
    ind = ind+100;
    for j = 1:100
        tof_data(ind(1,j),j) = 0;
    end
    imagesc(1:100,1:100,tof_data(101:end,:));
    xlabel('receiver')
    ylabel('transimitter')
    axis tight
    axis equal
    titlename = sprintf('Group1ERate%dIRate02Num3',i);
    axes;
    title(titlename);
    axis off;
    cd('H:\result\2018_04_27_analyzeLipidModel')
    myfilename = sprintf('2018_05_01_Group1ERate%dIRate02Num3_Compare-2',i);
    exportfig(myfilename,'png',[400,180]);
end

for i = 1:5
    figure;
    cd('H:\data\kwave\medium\2018_04_26_randomScatter')
    myfilename = sprintf('Group1ERate5IRate02Num%d-2', i);
    load(myfilename);
    subplot(1,2,1);
    imagesc(kgrid.x_vec*1000,kgrid.y_vec*1000,medium.sound_speed);
    axis equal
    axis tight
    hold on
    scatter(t_pos(2,:)*1000,t_pos(1,:)*1000,'rs','filled');
    xlabel('y-axis[mm]')
    ylabel('x-axis[mm]')
    hold off
    subplot(1,2,2);
    cd('H:\result\2018_04_27_analyzeLipidModel')
    myfilename = sprintf('2018_04_27_TOFdata_Group1ERate5IRate02Num%d-2',i);
    load(myfilename)
    [Min,ind]= min(tof_data(101:end,:));
    ind = ind+100;
    for j = 1:100
        tof_data(ind(1,j),j) = 0;
    end
    imagesc(1:100,1:100,tof_data(101:end,:));
    xlabel('receiver')
    ylabel('transimitter')
    axis tight
    axis equal
    titlename = sprintf('Group1ERate5IRate02Num%d',i);
    axes;
    title(titlename);
    axis off;
    cd('H:\result\2018_04_27_analyzeLipidModel')
    myfilename = sprintf('2018_05_01_Group1ERate5IRate02Num%d_Compare-2',i);
    exportfig(myfilename,'png',[400,180]);
end

for i = 1:5
    figure;
    cd('H:\data\kwave\medium\2018_04_26_randomScatter')
    myfilename = sprintf('Group1ERate5IRate%02dNum3-2', 2*i);
    load(myfilename);
    subplot(1,2,1);
    imagesc(kgrid.x_vec*1000,kgrid.y_vec*1000,medium.sound_speed);
    axis equal
    axis tight
    hold on
    scatter(t_pos(2,:)*1000,t_pos(1,:)*1000,'rs','filled');
    xlabel('y-axis[mm]')
    ylabel('x-axis[mm]')
    hold off
    subplot(1,2,2);
    cd('H:\result\2018_04_27_analyzeLipidModel')
    myfilename = sprintf('2018_04_27_TOFdata_Group1ERate5IRate%02dNum3-2',i*2);
    load(myfilename)
    [Min,ind]= min(tof_data(101:end,:));
    ind = ind+100;
    for j = 1:100
        tof_data(ind(1,j),j) = 0;
    end
    imagesc(1:100,1:100,tof_data(101:end,:));
    xlabel('receiver')
    ylabel('transimitter')
    axis tight
    axis equal
    titlename = sprintf('Group1ERate5IRate%02dNum3',i*2);
    axes;
    title(titlename);
    axis off;
    cd('H:\result\2018_04_27_analyzeLipidModel')
    myfilename = sprintf('2018_05_01_Group1ERate5IRate%02dNum3_Compare-2',i);
    exportfig(myfilename,'png',[400,180]);
end

for i = 2:2:20
    figure;
    cd('H:\data\kwave\medium\2018_04_26_randomScatter')
    myfilename = sprintf('Group1ERate%dIRate02Num3-2', i);
    load(myfilename);
    subplot(1,2,1);
    imagesc(kgrid.x_vec*1000,kgrid.y_vec*1000,medium.sound_speed);
    axis equal
    axis tight
    hold on
    scatter(t_pos(2,:)*1000,t_pos(1,:)*1000,'rs','filled');
    xlabel('y-axis[mm]')
    ylabel('x-axis[mm]')
    hold off
    subplot(1,2,2);
    cd('H:\result\2018_04_27_analyzeLipidModel')
    myfilename = sprintf('2018_04_27_TOFdata_Group1ERate%dIRate02Num3-2',i);
    load(myfilename)
    [Min,ind]= min(tof_data(101:end,:));
    ind = ind+100;
    for j = 1:100
        tof_data(ind(1,j),j) = 0;
    end
    imagesc(1:100,1:100,tof_data(101:end,:));
    xlabel('receiver')
    ylabel('transimitter')
    axis tight
    axis equal
    titlename = sprintf('Group1ERate%dIRate02Num3',i);
    axes;
    title(titlename);
    axis off;
    cd('H:\result\2018_04_27_analyzeLipidModel')
    myfilename = sprintf('2018_05_01_Group1ERate%dIRate02Num3_Compare-2',i);
    exportfig(myfilename,'png',[400,180]);
end
