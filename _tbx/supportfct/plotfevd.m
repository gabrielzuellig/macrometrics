function plotfevd(vd, printvars, path, printshocks, poslegend)

H = size(vd,1)-1;
n = size(vd,3);
nrow = ceil(n/3);
if isempty(poslegend)
    poslegend = 'NorthEast';
end
if isempty(printshocks)
    printshocks = {};
    for s = 1:size(vd,2)
        printshocks = [printshocks, ['shock', num2str(s)]];
    end
end

figure()
for vv = 1:n
    subplot(nrow, 3, vv)
    area(0:H, vd(:,:,vv));
    axis('tight')
    ylim([0 1])
    title(printvars(vv),'FontSize',14)
    if vv == 3 
        legend(printshocks, 'Location',poslegend,'FontSize',10)
        legend boxoff
    end
    set(gca,'GridLineStyle','-','Layer','bottom')
    if H <= 30
        set(gca,'XTickLabelMode', 'manual','XTickLabel',{'0','4','8','12','16','20','24'},'XTick',[0:4:24])
    else
        set(gca,'XTickLabelMode', 'manual','XTickLabel',{'0','12','24','36','48'},'XTick',[0:12:48])
    end
    set(gcf,'paperpositionmode','auto')
end
set(gcf, 'position', [0 0 800 nrow*200]);
print(gcf,'-depsc2','-loose',strcat(path,'fevd'));

