[sA_ind,sB_ind] = meshgrid(1:n_film,1:n_film);


subplot(2,1,1)


for ii = 1:length(sA_ind(:))
    plot_sticker_state(gr,sA_ind(ii),sB_ind(ii))
end
