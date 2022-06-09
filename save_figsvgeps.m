function save_figsvg(fn,saveit)

if saveit
f1=gcf()
savefig(f1,[fn '.fig'])
saveas(f1,[fn '.svg'],'svg')
saveas(f1,[fn '.eps'],'epsc')
end

end