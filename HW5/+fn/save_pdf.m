function save_pdf(h, name) 

% save as cropped pdf 
    set(h,'Units','Inches');
    pos = get(h,'Position');
    set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(h,name,'-dpdf','-r0')
    
end 