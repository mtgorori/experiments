tic
%%%% 初期パラメータ %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CT形状
ring_cx = 0; ring_cy = 0;% 中心位置
ring_radius = 50.e-3;% 半径
rr =  ring_radius;
%センサ設置
t_num = 16;%トランスデューサ数
t = linspace(0, ((t_num-1)/t_num)*2*pi, t_num);%センサ位置角
t_pos = zeros(2, t_num);%センサ位置
t_pos(1,:) = rr*cos(t)+ring_cx;
t_pos(2,:) = rr*sin(t)+ring_cy;
%セル設定
cell_num = 16;
cell_size = ring_radius*2 / (cell_num-1);
x_grid = -ring_radius: cell_size : ring_radius;
y_grid = x_grid;

p = zeros(2,length(x_grid)+length(y_grid));
C = zeros(t_num^2, cell_num^2);

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
            rm_ind = sqrt(p(1,:).^2 + p(2,:).^2)>rr;
            p(:,rm_ind) = [];
            
            %初期座標設定(xの値が最小な点）
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
                C(t_num*(ii-1)+jj,cell_num*(x_ind-1)+y_ind) = norm(p_init - p_neighbor);
                p_init = p_neighbor;
            end
        end
    end
end
A = sum(C,1);
figure;
plot(A);
axis tight
colorbar;
c = colorbar;
ylim([0 0.2])
hcbl = xlabel(c,'[m]');
xlabel('画素')
ylabel('経路密度')
exportfig( 'H:\result\matrix_rank_visualization\2018_03_02_matrix_rank_visualization_ring_16_1','png',[300,150]);