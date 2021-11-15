function plot_sticker_state(gr,sA_ind,sB_ind)
%PLOT_STICKER_STATE Summary of this function goes here
%   Detailed explanation goes here

X1 = [gr.A(sA_ind,1),gr.B(sB_ind,1)];
X2 = [gr.Dx(sA_ind,sB_ind),gr.Dx(sA_ind,sB_ind)];
Y1 = [gr.A(sA_ind,2),gr.B(sB_ind,2)];
Y2 = [gr.Dy(sA_ind,sB_ind),gr.Dy(sA_ind,sB_ind)];
plot([X1; X2], [Y1; Y2]);
end

