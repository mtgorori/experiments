tic
%%%% �����p�����[�^ %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%�Z���T�ݒu
t_size = 100.e-3;
t_num = 16;%�g�����X�f���[�T��
t_pos = zeros(2, t_num);%�Z���T�ʒu
t_pos(1,1:t_num/2) = -t_size/2:t_size/(t_num/2-1):t_size/2 ;%�f�q������������[m]
t_pos(2,1:t_num/2) = t_size/2;
t_pos(1,t_num/2+1:t_num) = t_pos(1,1:t_num/2);
t_pos(2,t_num/2+1:t_num) = -t_size/2;
%�Z���ݒ�
cell_num = 16;
p_size = 50.e-3;
cell_size = p_size*2 / (cell_num-1);
x_grid = -p_size: cell_size : p_size;
y_grid = x_grid;
p = zeros(2,length(x_grid)+length(y_grid));
C = zeros(cell_num, cell_num);

for ii = 1:t_num
    %���M�f�q�̍��W�ݒ�(x_tr,y_tr)
    x_tr = t_pos(1,ii);
    y_tr = t_pos(2,ii);
    pos_tr = [x_tr, y_tr];
    for jj = 1:t_num
        if (ii == jj) || (1<=ii)&&(ii<=t_num/2) && (1<=jj)&&(jj<=t_num/2) || (t_num/2<ii && t_num/2<jj)
            continue
        else
            %��M�f�q�̍��W�ݒ�(x_re,y_re)
            p_length = zeros(length(x_grid),length(y_grid));
            x_re = t_pos(1,jj);
            y_re = t_pos(2,jj);
            pos_re = [x_re, y_re];
            all_length = norm(pos_tr - pos_re);
            zz = zeros(length(x_grid),length(y_grid));
            %�f�q�����钼���̕�����
            x_line = ((x_tr-x_re) / (y_tr-y_re)) * (y_grid-y_re) + x_re;
            y_line = ((y_tr-y_re) / (x_tr-x_re)) * (x_grid-x_re) + y_re;
            %�e�f�q�ʒu���O���b�h�ɓ��Ă͂߂�
            [~,x_tr_index] = min(abs(x_tr - x_grid));
            [~,y_tr_index] = min(abs(y_tr - y_grid));
            [~,x_re_index] = min(abs(x_re - x_grid));
            [~,y_re_index] = min(abs(y_re - y_grid));
            %�ϕ��J�n�ʒu�C�I���ʒu�̌���i�����`�̈�):�o�H�����钷���`
            x_start = min(x_tr_index,x_re_index);
            x_end = max(x_tr_index,x_re_index);
            y_start = min(y_tr_index,y_re_index);
            y_end = max(y_tr_index,y_re_index);
            %�v�Z�̈撆�̊i�q�Ƃ̌�_���i�[(p)
            p(1,1:length(y_start:y_end)) = x_line(y_start:y_end);
            p(2,1:length(y_start:y_end)) = y_grid(y_start:y_end);
            p(1,length(y_start:y_end)+1:length(y_start:y_end)+length(x_start:x_end)) = x_grid(x_start:x_end);
            p(2,length(y_start:y_end)+1:length(y_start:y_end)+length(x_start:x_end)) = y_line(x_start:x_end);
            p(:,length(y_start:y_end)+length(x_start:x_end)+1:length(p)) = [];
            rm_ind = p(1,:)>t_size/2 | p(1,:)<-t_size/2 | p(2,:)>t_size/2 | p(2,:)<-t_size/2;
            p(:,rm_ind) = [];
            p=rmmissing(p,2);%NaN���o�����߁C�Ώ��D
            
            %�������W�ݒ�(x�̒l���ŏ��ȓ_�j
            [~,b] = min(p(1,:));
            p_init = p(:,b);
            while ~isempty(p)
                %������������Z���̌��o
                x_ind = find(abs(p_init(1,1) - x_grid) == 0, 1);
                y_ind = find(abs(p_init(2,1) - y_grid) == 0, 1);
                if isempty(x_ind)
                    [~,x_ind] = min(abs(p_init(1,1) - x_grid - cell_size/2));
                end
                if isempty(y_ind)
                    [~,y_ind] = min(abs(p_init(2,1) - y_grid - cell_size/2));
                end
                %�����̂�������̒[�_�̌��o
                p(:,b) = [];%�[�_�̏���
                [~,b] = min(p(1,:));%��������̒[�_�̃C���f�b�N�X���o
                p_neighbor = p(:,b);%��������̒[�_���o
                %��^�����Z�o
                C(y_ind,x_ind) = C(y_ind,x_ind)+ 1;
                p_init = p_neighbor;
            end
        end
    end
end
figure;
imagesc(x_grid*1000,y_grid*1000,C);
axis tight
axis equal
colorbar;
caxis([0 45]);
c = colorbar;
hcbl = xlabel(c,'�o�H���x');
xlabel('x����[mm]')
ylabel('y����[mm]')
exportfig( 'H:\result\matrix_rank_visualization\2018_03_02_matrix_rank_visualization_2board_16_2dim_2','png',[200,120]);