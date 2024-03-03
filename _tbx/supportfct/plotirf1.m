function plotirf1(irf, irfbands, printvars, path)

H = size(irf,1)-1;
n = size(irf,2);
lo = irfbands(:,:,1);
up = irfbands(:,:,2);
nrow = ceil(n/3);


figure()
for vv = 1:n
    maxim = max([irf(:,vv); lo(:,vv); up(:,vv)]);
    maxim = maxim + 0.1*abs(maxim);
    minim = min([irf(:,vv); lo(:,vv); up(:,vv)]);
    minim = minim - 0.1*abs(minim);
    subplot(nrow, 3, vv)
    plot(0:H, irf(:,vv), 'k', 'LineWidth', 2)
    hold on
    plot(0:H, zeros(H+1,1), 'k', 'LineWidth', 0.5)
    plot(0:H, lo(:,vv), '--k')
    plot(0:H, up(:,vv), '--k')
    grid on
    axis('tight')
    ylim([minim maxim])
    title(printvars(vv),'FontSize',14)
    set(gca,'GridLineStyle','-','Layer','bottom')
    if H <= 30
        set(gca,'XTickLabelMode', 'manual','XTickLabel',{'0','4','8','12','16','20','24'},'XTick',[0:4:24])
    else
        set(gca,'XTickLabelMode', 'manual','XTickLabel',{'0','12','24','36','48'},'XTick',[0:12:48])
    end
    set(gcf,'paperpositionmode','auto')
end
set(gcf, 'position', [0 0 800 nrow*200]);
print(gcf,'-dpng','-loose',strcat(path,'_irf.png'));
