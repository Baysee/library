function save_fig_png_svg_eps(fn)
% quick test to see if I can push back...
f1=gcf()
savefig(f1,[fn '.fig'])
saveas(f1,[fn '.png'],'png')
saveas(f1,[fn '.svg'],'svg')
saveas(f1,[fn '.eps'],'epsc')

end