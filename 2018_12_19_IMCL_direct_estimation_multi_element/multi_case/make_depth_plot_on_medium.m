%%%%%%%%%%
%correlation={num_assumed_SOS,num_assumed_depth,num_lateral}
%%%%%%%%%%

clear


case_name = [1,2,17,52];
for lll = 1:length(case_name)
    loadfilename = sprintf('H:/result/2018_12_19_IMCL_direct_estimation_multi_element/multi_case/2018_12_28_variousF&lateral/case%d2018_12_28_all_result',case_name(lll));
    load(loadfilename);
    for mm = 1:num_IMCL
        loadpath = sprintf('H:/data/kwave/result/2018_12_28_other_case_variousIMCL/case%d_IMCL%0.1f',case_name(lll),IMCL_rate(mm));
        load([loadpath,'/medium.mat'])
        cat_dst_path = sprintf('IMCL%d%%',IMCL_rate(mm));
        dst_path1 = [dst_path,cat_dst_path];
        load([dst_path1,'/result']);
        
        figure;
        imagesc(kgrid.x_vec*1000,kgrid.y_vec*1000,medium.sound_speed);
        hold on
        scatter(t_pos(2,1:t_num/2)*1000,t_pos(1,1:t_num/2)*1000,'rs','filled');
        for nn = 1:num_lateral
            if nn == ind_estimate_l
                continue
            end
            [~,ind_estimate_d] = max(max(correlation(:,:,nn)));
            scatter(assumed_point(2,ind_estimate_d,nn)*1000,assumed_point(1,ind_estimate_d,nn)*1000,'rs')
        end
        [~,ind_estimate_l] = max(max(max(correlation)));
        [~,ind_estimate_d] = max(max(correlation(:,:,ind_estimate_l)));
        scatter(assumed_point(2,ind_estimate_d,ind_estimate_l)*1000,...
            assumed_point(1,ind_estimate_d,ind_estimate_l)*1000,'rp')
        c = colorbar;
        c.Label.String = '[m/s]';
        set(gca,'YDir','normal');
        
        axis equal
        axis tight
        xlabel('x-axis[mm]')
        ylabel('y-axis[mm]')
        caxis([1450 1580]);
        xlim([5 20.7])
        ylim([-4.5 4.5])
        titlename = sprintf('IMCL= %d %%',IMCL_rate(mm));
        title(titlename)
        
        dst_path4 = [dst_path,'/estimated_depth'];
        if ~exist(dst_path4, 'dir')
            mkdir(dst_path4);
        end
        savefilename = sprintf('/IMCL%d%%',IMCL_rate(mm));
        savefig([dst_path4,savefilename,'.fig'])
        exportfig([dst_path4,savefilename],'png',[300,250])
        close gcf
        
    end
end
