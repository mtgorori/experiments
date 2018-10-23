%作成：2018/09/26-----------------------------------------------------------------------------------
%最短経路長の各素子ペアにおける個々の経路に関して，各媒質の専有面積比から算出した
%平均音速を配列v(送信素子番号,媒質番号,周波数)に格納するスクリプト．
%------------------------------------------------------------------------------------------------------
clear;close all
load("H:/data/kwave/result/2018_09_08_variousFrequency/statistics.mat",'tof_cell')
load("H:/data/kwave/param/param_2board.mat")
%% 初期設定
t_size = 40.e-3;
t_num = 200;%トランスデューサ数
t_pos = zeros(2, t_num);%センサ位置
t_pos(1,1:t_num/2) = -t_size/2:t_size/(t_num/2-1):t_size/2 ;%素子水平方向距離[m]
t_pos(2,1:t_num/2) = t_size/2;
t_pos(1,t_num/2+1:t_num) = t_pos(1,1:t_num/2);
t_pos(2,t_num/2+1:t_num) = -t_size/2;
[r_num,t_num,model_num,freq_num] = size(tof_cell);
leng = zeros(r_num,t_num);
for i = 1:r_num
    for j = 1:t_num
        leng(i,j) = sqrt((t_pos(1,i)-t_pos(1,j))^2+(t_pos(2,i)-t_pos(2,j))^2);
    end
end
%% lt(=v)データの作成
v = zeros(t_num,model_num,freq_num);
for ii = 1:freq_num
    for jj = 1:model_num
        for kk = 1:t_num
            v(kk,jj,ii) = leng(100+kk,kk)/tof_cell(100+kk,kk,jj,ii);
        end
    end
end
figure;
plot(v(:,1,2));
figure;
waterfall(v(:,:,1)');
xlabel('素子番号')
ylabel('媒質番号')
zlabel('平均音速[m/s]')
exportfig("H:\result\2018_09_26_ltplot\lt_from_signal_waterfall",'png',[400,400])
figure;
imagesc(v(:,:,1)');
xlabel('素子番号')
ylabel('媒質番号')
axis tight
axis square
c = colorbar;
c.Label.String = '平均音速[m/s]';
exportfig("H:\result\2018_09_26_ltplot\lt_from_signal_imagesc",'png',[300,300])
%% 参照用のlt(=v)データの作成
v_reference = zeros(t_num,model_num);
% 適切な参照範囲を探す．
displace_t_m = abs(t_pos(1,2) - t_pos(1,1));%[m]
displace_t = displace_t_m / param.grid.dx;%[t_num]
% 一素子につきグリッド2個分の画素値を平均化して参照用のlt(=v)データの作成を行う．
load("H:/data/kwave/medium/2018_08_10_realisticScatter/case1.mat")
figure;
imagesc(kgrid.x_vec*1000,kgrid.y_vec*1000,medium.sound_speed);
ax = gca;
ax.YDir = 'normal';
xlabel('x軸[mm]')
ylabel('y軸[mm]')
exportfig("H:\result\2018_09_26_ltplot\medium_case1",'png',[300,300]);
%  素子がどのグリッドに位置しているかを求める．
ind_t_pos = zeros(2,200);
[~,ind_t_pos(1,:)] = min(abs(t_pos(1,:)-kgrid.y_vec));
[~,ind_t_pos(2,:)] = min(abs(t_pos(2,:)-kgrid.x_vec));
ref_medium = medium;
%どのグリッドを参照して参照用のltデータを作成しているかを確認する．また，試行的に一番目のモデルに関して参照用のデータを作成する．
for i = 1:t_num
    v_reference(i,1) = mean2(ref_medium.sound_speed(ind_t_pos(1,i)-1:ind_t_pos(1,i)+1,ind_t_pos(2,100+i):ind_t_pos(2,i)));
    ref_medium.sound_speed(ind_t_pos(1,i)-1:ind_t_pos(1,i)+1,ind_t_pos(2,100+i):ind_t_pos(2,i)) = 0;
end
figure;
imagesc(kgrid.x_vec*1000,kgrid.y_vec*1000,ref_medium.sound_speed);
ax = gca;
ax.YDir = 'normal';
xlabel('x軸[mm]')
ylabel('y軸[mm]')
exportfig("H:\result\2018_09_26_ltplot\refer_range_medium_case1",'png',[300,300]);
v_reference = zeros(t_num,model_num);
% 画素値平均化を行う
for jj = 1:model_num
    loadfilename = sprintf("H:/data/kwave/medium/2018_08_10_realisticScatter/case%d.mat",jj);
    load(loadfilename);
    ref_medium = medium;
    for ii = 1:t_num
        v_reference(ii,jj) = mean2(ref_medium.sound_speed(ind_t_pos(1,ii)-1:ind_t_pos(1,ii)+1,ind_t_pos(2,100+ii):ind_t_pos(2,ii)));
    end
end
v_reference = repmat(v_reference,1,1,6);
figure;
waterfall(v_reference(:,:,1)');
xlabel('素子番号')
ylabel('媒質番号')
exportfig("H:\result\2018_09_26_ltplot\lt_from_medium_waterfall",'png',[300,300]);
figure;
imagesc(v_reference(:,:,1)');
xlabel('素子番号')
ylabel('媒質番号')
axis tight
axis square
c = colorbar;
c.Label.String = '平均音速[m/s]';
exportfig("H:\result\2018_09_26_ltplot\lt_from_medium_imagesc",'png',[300,300])
%% vプロットの比較
%water-fallプロット
figure;
subplot(1,2,1);
waterfall(v_reference(:,:,1)');
xlabel('素子番号')
ylabel('媒質番号')
zlabel('平均音速[m/s]')
title({'素子ペアごとの平均音速値';'(媒質データから作成)'})
axis square
subplot(1,2,2);
waterfall(v(:,:,1)');
xlabel('素子番号')
ylabel('媒質番号')
zlabel('平均音速[m/s]')
title({'素子ペアごとの平均音速値';'(RFデータから作成)'})
axis square
exportfig("H:\result\2018_09_26_ltplot\lt_compared_waterfall",'png',[800,400])
%imagescプロット
figure;
subplot(1,2,1);
imagesc(v_reference(:,:,1)');
xlabel('素子番号')
ylabel('媒質番号')
title({'素子ペアごとの平均音速値';'(媒質データから作成)'})
axis square
colorbar
caxis([1450 1580])
ax = gca;
ax.YDir = 'normal';
subplot(1,2,2);
imagesc(v(:,:,1)');
xlabel('素子番号')
ylabel('媒質番号')
title({'素子ペアごとの平均音速値';'(RFデータから作成)'})
axis square
colorbar
caxis([1450 1580])
ax = gca;
ax.YDir = 'normal';
exportfig("H:\result\2018_09_26_ltplot\lt_compared_imagesc",'png',[800,400])
%% 各種データ保存