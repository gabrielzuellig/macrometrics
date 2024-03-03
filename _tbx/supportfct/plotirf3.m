function plotirf3(irf1, irfbands1, irf2, irfbands2, irf3, irfbands3, printvars, path, printstate, poslegend)

H = size(irf1,1)-1;
n = size(irf1,2);
lo1 = irfbands1(:,:,1);
up1 = irfbands1(:,:,2);
lo2 = irfbands2(:,:,1);
up2 = irfbands2(:,:,2);
lo3 = irfbands3(:,:,1);
up3 = irfbands3(:,:,2);
nrow = ceil(n/3);
if isempty(poslegend)
    poslegend = 'NorthEast';
end


figure()
for vv = 1:n
    maxim = max([irf1(:,vv); irf2(:,vv); irf3(:,vv); ...
                 up1(:,vv); up2(:,vv); up3(:,vv)]);
    maxim = maxim + 0.1*abs(maxim);
    minim = min([irf1(:,vv); irf2(:,vv); irf3(:,vv); ...
                 lo1(:,vv); lo2(:,vv); lo3(:,vv)]);
    minim = minim - 0.1*abs(minim);
    subplot(nrow, 3, vv)
    shadedplot(0:H,lo1(:,vv)',up1(:,vv)',[.9 .9 .9],[.9 .9 .9] );
    hold on
    p1=plot(0:H, irf1(:,vv), 'k', 'LineWidth', 2);
    plot(0:H, zeros(H+1,1), 'k', 'LineWidth', 0.5)
    p2=plot(0:H, irf2(:,vv), '--b', 'LineWidth', 2);
    plot(0:H, lo2(:,vv), '--b')
    plot(0:H, up2(:,vv), '--b')
    p3=plot(0:H, irf3(:,vv), ':r', 'LineWidth', 2);
    plot(0:H, lo3(:,vv), ':r')
    plot(0:H, up3(:,vv), ':r')
    grid on
    axis('tight')
    ylim([minim maxim])
    title(printvars(vv),'FontSize',14)
    if vv == 3 && ~isempty(printstate)
        legend([p1 p2 p3], printstate, 'Location',poslegend,'FontSize',10)
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
print(gcf,'-dpng','-loose',strcat(path,'_irf.png'));