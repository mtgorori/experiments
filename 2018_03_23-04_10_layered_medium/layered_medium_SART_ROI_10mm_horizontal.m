tic
%make_layered_medium
%概要説明：筋肉，脂肪，筋肉の三層構造(トランスデューサに対して水平向き）．脂肪層の厚さが線形に増加すると，ROI中の平均音速も線形に増加する．
%目的：このデータを基準にして音速分布再構成や音速推定を行って，真値との差を評価する．
%   詳細説明をここに記述
%　input:Layer.thicknessー脂肪層の厚さ[mm]
%　output:Medium.sound_speedー媒質の音速分布[m/s]
%               Medium.densityー媒質の密度分布[kg/m3]
%               Layer.sound_speed_aveーROI中の平均音速[m/s]
% 音速
% v_muscle = 1580;%筋肉の音速[m/s]
% v_fat = 1450;%脂肪の音速[m/s]
% 密度
% den_muscle = 1040;%筋肉の密度[kg/m3]
% den_fat = 920;%脂肪の密度[kg/m3]
clear 
close all
% 初期パラメータ
Sensor.sizeTotal = 100.e-3;%トランスデューサ長[m]
Sensor.num = 256;%トランスデューサ数[-]
thickness_Min = 2;%[mm]
thickness_Max = 10;%[mm]
iteration = 30;
% グリッド設定
[ Grid ] = make_grid( Sensor );
% センサ設定
[ Sensor ] = make_sensor_2board( Sensor );

us = zeros(length(Grid.x),length(Grid.y),iteration,thickness_Max);
RMSE = zeros(1,iteration,thickness_Max);

for thickness = thickness_Min:2:thickness_Max%脂肪層の厚さ[mm]
    % 媒質設定
    [ Medium, Layer ] = make_layered_medium_horizontal( thickness, Grid );
    % 投影データ作成
    [ projection ] = getProjectionData( Sensor, Grid, Medium );
    % 音速再構成
    pjt_est = zeros(Sensor.num,Sensor.num);
    p = zeros(2,length(Grid.x)+length(Grid.y));%交点座標群
    reI = zeros(length(Grid.x),length(Grid.y));
    reI_store = zeros(length(Grid.x),length(Grid.y),iteration);
    det_I = zeros(length(Grid.x),length(Grid.y));
    for it = 1:iteration
        count = ones(length(Grid.x));
        for ii = 1:Sensor.num
            %送信素子の座標設定(x_tr,y_tr)
            x_tr = Sensor.pos(1,ii);
            y_tr = Sensor.pos(2,ii);
            pos_tr = [x_tr, y_tr];
            for jj = 1:Sensor.num
                if (ii == jj) || (1<=ii)&&(ii<=Sensor.num/2) && (1<=jj)&&(jj<=Sensor.num/2) || (Sensor.num/2<ii && Sensor.num/2<jj)
                    continue
                else
                    %受信素子の座標設定(x_re,y_re)
                    p_length = zeros(length(Grid.x),length(Grid.y));
                    x_re = Sensor.pos(1,jj);
                    y_re = Sensor.pos(2,jj);
                    pos_re = [x_re, y_re];
                    all_length = norm(pos_tr - pos_re);
                    zz = zeros(length(Grid.x),length(Grid.y));
                    %素子が張る直線の方程式
                    x_line = ((x_tr-x_re) / (y_tr-y_re)) * (Grid.y-y_re) + x_re;
                    y_line = ((y_tr-y_re) / (x_tr-x_re)) * (Grid.x-x_re) + y_re;
                    %各素子位置をグリッドに当てはめる
                    [~,x_tr_index] = min(abs(x_tr - Grid.x));
                    [~,y_tr_index] = min(abs(y_tr - Grid.y));
                    [~,x_re_index] = min(abs(x_re - Grid.x));
                    [~,y_re_index] = min(abs(y_re - Grid.y));
                    %積分開始位置，終了位置の決定（長方形領域):経路が張る長方形
                    x_start = min(x_tr_index,x_re_index);
                    x_end = max(x_tr_index,x_re_index);
                    y_start = min(y_tr_index,y_re_index);
                    y_end = max(y_tr_index,y_re_index);
                    %計算領域中の格子との交点を格納(p)
                    p(1,1:length(y_start:y_end)) = x_line(y_start:y_end);
                    p(2,1:length(y_start:y_end)) = Grid.y(y_start:y_end);
                    p(1,length(y_start:y_end)+1:length(y_start:y_end)+length(x_start:x_end)) = Grid.x(x_start:x_end);
                    p(2,length(y_start:y_end)+1:length(y_start:y_end)+length(x_start:x_end)) = y_line(x_start:x_end);
                    p(:,length(y_start:y_end)+length(x_start:x_end)+1:length(p)) = [];
                    rm_ind = p(1,:)>Sensor.sizeTotal/2 | p(1,:)<-Sensor.sizeTotal/2 | p(2,:)>Sensor.sizeTotal/2 | p(2,:)<-Sensor.sizeTotal/2;
                    p(:,rm_ind) = [];
                    p=rmmissing(p,2);%NaNが出たため，対処．
                    %初期座標設定(xの値が最小な点）
                    if max(p(1,:)) - min(p(1,:)) < Grid.size % argx(min(p))が正しく判定できないことを防ぐため．
                        [~,b] = min(p(2,:));
                        p_init = p(:,b);
                        while ~isempty(p)
                            %線分が属するセルの検出
                            x_ind = find(abs(p_init(1,1) - Grid.x) == 0, 1);
                            y_ind = find(abs(p_init(2,1) - Grid.y) == 0, 1);
                            if isempty(x_ind)
                                [~,x_ind] = min(abs(p_init(1,1) - Grid.x - Grid.size/2));
                            end
                            if isempty(y_ind)
                                [~,y_ind] = min(abs(p_init(2,1) - Grid.y - Grid.size/2));
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
                            x_ind = find(abs(p_init(1,1) - Grid.x) == 0, 1);
                            y_ind = find(abs(p_init(2,1) - Grid.y) == 0, 1);
                            if isempty(x_ind)
                                [~,x_ind] = min(abs(p_init(1,1) - Grid.x - Grid.size/2));
                            end
                            if isempty(y_ind)
                                [~,y_ind] = min(abs(p_init(2,1) - Grid.y - Grid.size/2));
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
                pjt_act = projection(jj,ii);
                count = count + zz;
                Det = ((pjt_act - pjt_est(jj,ii))/all_length).*p_length;
                det_I = det_I +Det;
            end
        end
        reI = reI + det_I./count;
        reI_store(:,:,it) = reI;
        us(:,:,it,thickness) = Grid.size./(reI+Grid.size/Medium.v0);
        det_I = zeros(length(Grid.x),length(Grid.y));
        RMSE(1,it,thickness) = sqrt(immse(us(:,:,it,thickness), Medium.sound_speed));
    end
end
toc