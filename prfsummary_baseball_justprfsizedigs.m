
% Same but for digits

close all
close all hidden

whichSeparation = DIGDEXbloop;

clear g
figure('Position',[100 100 1000 500])
g(1,1) = gramm('x',whichSeparation,'y',[D2;D3;D4;D5],'color',whichSeparation);
g(1,1).geom_jitter()
%g(1,1).stat_bin('geom','line')
g(1,1).no_legend()
g(1,1).set_title('Within digit direction')

g(1,2) = gramm('x',whichSeparation,'y',[D2y;D3y;D4y;D5y],'color',whichSeparation);
g(1,2).geom_jitter()
g(1,2).no_legend()
g(1,2).set_title('Between digit direction')



g.set_names('x',[],'y', 'pRF size')
g.set_text_options('Font','Helvetica', 'base_size', 14)
g.set_point_options('base_size',8)
g.axe_property('XGrid','on','YGrid','on')
if thisModel == 2
    g.axe_property('YLim',[0 4])
else
    g.axe_property('YLim',[0 5])
end
g.set_order_options('x',0)
g.draw

g(1,1).update
g(1,1).stat_boxplot('width', 2, 'dodge', 0.5, 'alpha', 0.5, 'linewidth', 2, 'drawoutlier',0)
g(1,1).no_legend()


g(1,2).update
g(1,2).stat_boxplot('width', 2, 'dodge', 0.5, 'alpha', 0.5, 'linewidth', 2, 'drawoutlier',0)
g(1,2).no_legend()

g.draw
filename = sprintf(['DIG_sue_compare_' sep '_' datafilenames{thisModel} '.pdf'],'%s%s');
g.export('file_name',filename, ...
    'export_path',...
    '/Users/ppzma/The University of Nottingham/Michael_Sue - pRF/pRF_2021/',...
    'file_type','pdf')





