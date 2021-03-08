function increase_ylim

    ylims  = get(gca, 'ylim');
    yd     = ylims(2) - ylims(1); 
    set(gca, 'ylim', [ylims(1) - 0.2*yd, ylims(2) + 0.2*yd  ]); 

end 