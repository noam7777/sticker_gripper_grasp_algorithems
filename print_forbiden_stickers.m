
[sA_ind,sB_ind] = find(gr.LA_loose_cond_geo);





for ii = 1:length(sA_ind)
    plot_sticker_state(gr,sA_ind(ii),sB_ind(ii))
end
