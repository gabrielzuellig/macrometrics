function plothd(hd, time, printvars, path, shocknames, poslegend)

t = size(hd,1);
n = size(hd,3);
if isempty(poslegend)
    poslegend = 'NorthEast';
end
if isempty(shocknames)
    shocknames = {};
    for s = 1:size(vd,2)
        shocknames = [shocknames, ['shock', num2str(s)]];
    end
end


for vv = 1:n
    figure()
    X = squeeze(hd(:,:,vv));
    H(1,:) = area(time, (X).*(X>0),'LineStyle','none'); 
    hold on
    H(2,:) = area(time, (X).*(X<0),'LineStyle','none'); 
    for ii=1:size(H,2)
        H(2,ii).FaceColor = H(1,ii).FaceColor;
        H(2,ii).EdgeColor = H(2,ii).FaceColor;
    end
    h = plot(time, sum(squeeze(hd(:,:,vv)),2),'-k','LineWidth',2);
    grid on
    axis('tight')
    limits = ylim;
    minim = limits(1) - 0.1*abs(limits(1));
    maxim = limits(2) + 0.1*abs(limits(2));
    if vv == 1
        minim = limits(1) - 0.2*abs(limits(1));
        maxim = limits(2) + 0.2*abs(limits(2));
        legend(shocknames, 'Location',poslegend,'FontSize',10)
        legend boxoff
    end
    ylim([minim maxim])
    title(printvars(vv),'FontSize',14)
    set(gca,'GridLineStyle','-','Layer','bottom')
    set(gcf,'paperpositionmode','auto')
    set(gcf, 'position', [0 0 800 200]);
    print(gcf,'-depsc2','-loose',strcat(path,'varhd_', num2str(vv)));
end

