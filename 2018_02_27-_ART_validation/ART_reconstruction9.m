%%%%TH metod%%%%%%%%%%

clear
close all
tic
%%%% 初期パラメータ %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CT形状
ring_cx = 0; ring_cy = 0;% 中心位置
ring_radius = 50.e-3;% 半径
rr =  ring_radius;
%音速
v_water = 1540;%水の音速[m/s]
v_fat = 1420;%脂肪の音速[m/s]
%脂肪領域形状
fat_cx = 20.e-3; fat_cy = 0; % 中心位置
fat_radius = 10.e-3;
fd = fat_radius;
%センサ設置
t_num = 256;%トランスデューサ数
t = linspace(0, ((t_num-1)/t_num)*2*pi, t_num);%センサ位置角
t_pos = zeros(2, t_num);%センサ位置
t_pos(1,:) = rr*cos(t)+ring_cx;
t_pos(2,:) = rr*sin(t)+ring_cy;
%セル設定
cell_num = 256;
cell_size = ring_radius*2 / cell_num;
x_grid = -ring_radius: cell_size : ring_radius;
y_grid = x_grid;
%正解音速分布
v_dist = v_water*ones(length(x_grid),length(y_grid));
for i = 1:length(x_grid)
    for j = 1:length(y_grid)
        if (x_grid(i)-fat_cx)^2 + (y_grid(j)-fat_cy)^2 <= fd^2
            v_dist(j,i) = v_fat;
        end
    end
end
%%%% 音速分布 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% 実測投影データ呼び出し%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('../result/tof_aic_win_200.mat','tof_map')
%%%% 音速再構成 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
iteration = 10;
pjt_est = zeros(t_num,t_num);
p = zeros(2,length(x_grid)+length(y_grid));
reI = zeros(length(x_grid),length(y_grid));%再構成画像
reI_store = zeros(length(x_grid),length(y_grid),iteration);
us = zeros(length(x_grid),length(y_grid),iteration);
MSE = zeros(1,iteration);

for it = 1:iteration
    for ii = 1:t_num
        %送信素子の座標設定(x_tr,y_tr)
        x_tr = t_pos(1,ii);
        y_tr = t_pos(2,ii);
        pos_tr = [x_tr, y_tr];
        for jj = 1:t_num
            if ii == jj
                continue
            else
                %受信素子の座標設定(x_re,y_re)
                p_length = zeros(length(x_grid),length(y_grid));
                x_re = t_pos(1,jj);
                y_re = t_pos(2,jj);
                pos_re = [x_re, y_re];
                all_length = norm(pos_tr - pos_re);
                zz = zeros(length(x_grid),length(y_grid));
                %素子が張る直線の方程式
                x_line = ((x_tr-x_re) / (y_tr-y_re)) * (y_grid-y_re) + x_re;
                y_line = ((y_tr-y_re) / (x_tr-x_re)) * (x_grid-x_re) + y_re;
                %各素子位置をグリッドに当てはめる
                [~,x_tr_index] = min(abs(x_tr - x_grid));
                [~,y_tr_index] = min(abs(y_tr - y_grid));
                [~,x_re_index] = min(abs(x_re - x_grid));
                [~,y_re_index] = min(abs(y_re - y_grid));
                %積分開始位置，終了位置の決定（長方形領域):経路が張る長方形
                x_start = min(x_tr_index,x_re_index);
                x_end = max(x_tr_index,x_re_index);
                y_start = min(y_tr_index,y_re_index);
                y_end = max(y_tr_index,y_re_index);
                %計算領域中の格子との交点を格納(p)
                p(1,1:length(y_start:y_end)) = x_line(y_start:y_end);
                p(2,1:length(y_start:y_end)) = y_grid(y_start:y_end);
                p(1,length(y_start:y_end)+1:length(y_start:y_end)+length(x_start:x_end)) = x_grid(x_start:x_end);
                p(2,length(y_start:y_end)+1:length(y_start:y_end)+length(x_start:x_end)) = y_line(x_start:x_end);
                p(:,length(y_start:y_end)+length(x_start:x_end)+1:length(p)) = [];
                rm_ind = find(sqrt(p(1,:).^2 + p(2,:).^2)>rr);
                p(:,rm_ind) = [];
                %初期座標設定(xの値が最小な点）
                if max(p(1,:)) - min(p(1,:)) < 3*cell_size % argx(min(p))が正しく判定できないことを防ぐため．
                    [~,b] = min(p(2,:));
                    p_init = p(:,b);
                    while ~isempty(p)
                        %線分が属するセルの検出
                        x_ind = find(abs(p_init(1,1) - x_grid) == 0, 1);
                        y_ind = find(abs(p_init(2,1) - y_grid) == 0, 1);
                        if isempty(x_ind)
                            [~,x_ind] = min(abs(p_init(1,1) - x_grid - cell_size/2));
                        end
                        if isempty(y_ind)
                            [~,y_ind] = min(abs(p_init(2,1) - y_grid - cell_size/2));
                        end
                        %線分のもう一方の端点の検出
                        p(:,b) = [];%端点の除去
                        [~,b] = min(p(2,:));%もう一方の端点のインデックス検出
                        p_neighbor = p(:,b);%もう一方の端点検出
                        %寄与音速算出
                        p_length(y_ind,x_ind) = norm(p_init - p_neighbor);
%                         pjt_est(ii,jj) = pjt_est(ii,jj) + p_length(x_ind,y_ind)*reI(ii,jj);
                        p_init = p_neighbor;
                        zz(y_ind,x_ind) = 1;
                    end
                else
                    [~,b] = min(p(1,:));
                    p_init = p(:,b);
                    while ~isempty(p)
                        %線分が属するセルの検出
                        x_ind = find(abs(p_init(1,1) - x_grid) == 0, 1);
                        y_ind = find(abs(p_init(2,1) - y_grid) == 0, 1);
                        if isempty(x_ind)
                            [~,x_ind] = min(abs(p_init(1,1) - x_grid - cell_size/2));
                        end
                        if isempty(y_ind)
                            [~,y_ind] = min(abs(p_init(2,1) - y_grid - cell_size/2));
                        end
                        %線分のもう一方の端点の検出
                        p(:,b) = [];%端点の除去
                        [~,b] = min(p(1,:));%もう一方の端点のインデックス検出
                        p_neighbor = p(:,b);%もう一方の端点検出
                        %寄与音速算出
                        p_length(y_ind,x_ind) = norm(p_init - p_neighbor);
                        p_init = p_neighbor;
                        zz(y_ind,x_ind) = 1;
                    end
                end
            end
            pjt_est(jj,ii) = sum(sum(reI.*zz));%(ii,jj)はiiが行，すなわちy方向，jjが列，すなわちx方向なので転置して(x,y)に置き換える必要があった．
            pjt_act = tof_map(jj,ii);
            det_I = (pjt_act - pjt_est(jj,ii))/all_length;
            reI = reI + det_I.*p_length;
        end
    end
    reI_store(:,:,it) = reI;
    us(:,:,it) = cell_size./(reI+cell_size/v_water);
    MSE(1,it) = immse(us(:,:,it), v_dist);
    if (it>1)&&((MSE(1,it-1)-MSE(1,it))<0.1)
          It = it-1;    
        break
    end
end
if It < iteration
    it = It;
end
toc
%%%% 画質評価 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
plot(MSE(1,1:it));
xlabel('更新回数')
ylabel('平均二乗誤差[m/s]')
axis tight
exportfig('../result/ART_reconstruction_AIC200_MSE','png',[250,200]);
%%%% プロファイル表示　%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
plot(x_grid*1000,us(t_num/2, :,it))
hold on
plot(x_grid*1000,v_dist(t_num/2,:),'r');
xlabel('x方向[mm]')
ylabel('音速[m/s]')
legend('再構成音速','正解音速','Location','southwest')
hold off
exportfig('../result/ART_reconstruction_AIC200_profile','png',[300,300]);
%%%% 描画用その1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
subplot(1,2,1);
imagesc(x_grid*1000, y_grid*1000, v_dist);
hold on
scatter(t_pos(1,:)*1000,t_pos(2,:)*1000,3,'r');
hold off
colorbar;
set(gca,'YDir','normal');
xlabel('x方向[mm]')
ylabel('y方向[mm]')
caxis([1350 1600])
axis equal;
axis tight
subplot(1,2,2)
imagesc(x_grid*1000,y_grid*1000,us(:,:,it));
colorbar;
axis equal;
axis tight;
set(gca,'YDir','normal');
xlabel('x方向[mm]')
ylabel('y方向[mm]')
caxis([1350 1600])
toc
exportfig('../result/ART_reconstruction_AIC200_com','png',[800,400]);
%%%% 描画用その2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%いったん収束し，その後発散する様子を示す．
figure
for i = 1:4
subplot(2,3,i);
imagesc(x_grid*1000,y_grid*1000,us(:,:,i));
c=colorbar;
c.Label.String = '[m/s]';
axis equal;
axis tight;
set(gca,'YDir','normal');
title(['iteration count:',num2str(i),''])
xlabel('x方向[mm]')
ylabel('y方向[mm]')
caxis([1350 1600])
end
i = 10;
subplot(2,3,5);
imagesc(x_grid*1000,y_grid*1000,us(:,:,i));
c=colorbar;
c.Label.String = '[m/s]';
axis equal;
axis tight;
set(gca,'YDir','normal');
title(['iteration count:',num2str(i),''])
xlabel('x方向[mm]')
ylabel('y方向[mm]')
caxis([1350 1600])
i = 15;
subplot(2,3,6);
imagesc(x_grid*1000,y_grid*1000,us(:,:,i));
c=colorbar;
c.Label.String = '[m/s]';
axis equal;
axis tight;
set(gca,'YDir','normal');
title(['iteration count:',num2str(i),''])
xlabel('x方向[mm]')
ylabel('y方向[mm]')
caxis([1350 1600])
exportfig("H:\result\2018_02_27-_ART_validation\ART_reconstruction5_iter",'png',[1000,450]);
% exportfig('../result/ART_reconstruction_AIC200_iter','png',[1000,450]);
% %%%% 描画用その3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%正解音速分布と再構成音速分布の差のヒストグラムを表示させる．
% figure
% histogram(v_dist - us(:,:,it-1),'Normalization','probability');
%%%%実際はこのような考え方で組んだが，転置処理に時間がかかる
%%%%（iteration:2, 303.320936秒）ので，それをp_length内などで行ってしまうというやりかたを上でやった．（iteration:2, 261.735875 秒)
% 
%  if max(p(1,:)) - min(p(1,:)) < 3*cell_size % argx(min(p))が正しく判定できないことを防ぐため．
%                     [~,b] = min(p(2,:));
%                     p_init = p(:,b);
%                     while ~isempty(p)
%                         %線分が属するセルの検出
%                         x_ind = find(abs(p_init(1,1) - x_grid) == 0, 1);
%                         y_ind = find(abs(p_init(2,1) - y_grid) == 0, 1);
%                         if isempty(x_ind)
%                             [~,x_ind] = min(abs(p_init(1,1) - x_grid - cell_size/2));
%                         end
%                         if isempty(y_ind)
%                             [~,y_ind] = min(abs(p_init(2,1) - y_grid - cell_size/2));
%                         end
%                         %線分のもう一方の端点の検出
%                         p(:,b) = [];%端点の除去
%                         [~,b] = min(p(2,:));%もう一方の端点のインデックス検出
%                         p_neighbor = p(:,b);%もう一方の端点検出
%                         %寄与音速算出
%                         p_length(x_ind,y_ind) = norm(p_init - p_neighbor);
% %                         pjt_est(ii,jj) = pjt_est(ii,jj) + p_length(x_ind,y_ind)*reI(ii,jj);
%                         p_init = p_neighbor;
%                         zz(x_ind,y_ind) = 1;
%                     end
%                 else
%                     [~,b] = min(p(1,:));
%                     p_init = p(:,b);
%                     while ~isempty(p)
%                         %線分が属するセルの検出
%                         x_ind = find(abs(p_init(1,1) - x_grid) == 0, 1);
%                         y_ind = find(abs(p_init(2,1) - y_grid) == 0, 1);
%                         if isempty(x_ind)
%                             [~,x_ind] = min(abs(p_init(1,1) - x_grid - cell_size/2));
%                         end
%                         if isempty(y_ind)
%                             [~,y_ind] = min(abs(p_init(2,1) - y_grid - cell_size/2));
%                         end
%                         %線分のもう一方の端点の検出
%                         p(:,b) = [];%端点の除去
%                         [~,b] = min(p(1,:));%もう一方の端点のインデックス検出
%                         p_neighbor = p(:,b);%もう一方の端点検出
%                         %寄与音速算出
%                         p_length(x_ind,y_ind) = norm(p_init - p_neighbor);
% %                         pjt_est(ii,jj) = pjt_est(ii,jj) + reI(x_ind,y_ind);
%                         p_init = p_neighbor;
%                         zz(x_ind,y_ind) = 1;
%                     end
%                 end
% pjt_est(jj,ii) = sum(sum(reI.*zz.'));%(ii,jj)はiiが行，すなわちy方向，jjが列，すなわちx方向なので転置して(x,y)に置き換える必要があった．
%             pjt_act = tof_map(jj,ii);
%             det_I = (pjt_act - pjt_est(jj,ii))/all_length;
%             reI = reI + det_I.*p_length.';
%         end