function save_figpng(fn)

f1=gcf()
savefig(f1,[fn '.fig'])
saveas(f1,[fn '.png'],'png')


end