%%%%%%%%%%
%correlation={num_assumed_SOS,num_assumed_depth,num_lateral}
%%%%%%%%%%

clear


case_name = [1,2,17,52];
for lll = 1:length(case_name)
    loadfilename = sprintf('H:/result/2018_12_19_IMCL_direct_estimation_multi_element/multi_case/2018_12_28_variousF&lateral/case%d2018_12_28_all_result',case_name(lll));
    load(loadfilename);
    for mm = 1:num_IMCL
        cat_dst_path = sprintf('IMCL%d%%',IMCL_rate(mm));
        dst_path1 = [dst_path,cat_dst_path];
        load([dst_path1,'/result']);
        
        figure;
        [~,ind_estimate_l] = max(max(max(correlation)));
        [~,ind_estimate_d] = max(max(correlation(:,:,ind_estimate_l)));
        subplot(1,2,1)
        imagesc(assumed_depth*1e3,assumed_IMCL_rate,correlation(:,:,ind_estimate_l));
        axis tight
        axis square
        hold on
        scatter(assumed_depth(ind_estimate_d)*1e3,assumed_IMCL_rate(ind_estimate_v),100,'rp')
        xlabel('x-axis[mm]')
        ylabel('estimated IMCL content [%]')
        title('correlation')
        hold off
        
        subplot(1,2,2)
        plot(correlation(:,ind_estimate_d,ind_estimate_l));
        ax2 = subplot(1,2,2);
        ylim manual
        ylim(ax2,[-0.4 0.4])
        axis tight
        axis square
        xlabel('estimated IMCL content [%]')
        ylabel('correlation')
        hold on
        scatter(ind_estimate_v,correlation(ind_estimate_v,ind_estimate_d,ind_estimate_l),100,'rp');
        hold off
        
        axes;
        titlename = sprintf('IMCL= %d %%',IMCL_rate(mm));
        title(titlename)
        axis off;
        dst_path4 = [dst_path,'/correlation_plot'];;
        if ~exist(dst_path4, 'dir')
            mkdir(dst_path4);
        end
        savefilename = sprintf('/IMCL%d%%',IMCL_rate(mm));
        savefig([dst_path4,savefilename,'.fig'])
        exportfig([dst_path4,savefilename],'png',[400,190])
        close gcf
        
    end
end
