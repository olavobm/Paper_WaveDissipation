function disp_figstart()
%% Print on the screen my usual code for starting a new figure.
% (then I can just copy and paste on the editor).


%
disp('%')
disp('hfig = figure;')
disp('hfig.Units = ''normalized'';')
disp('hfig.Position = [0.2, 0.2, 0.5, 0.5];')
disp('%')
disp('haxs = makeSubPlots(0.1, 0.1, 0.1, ...')
disp('                    0.1, 0.1, 0.1, 1, 1);')
disp('hold(haxs, ''on'')')
disp(' ')
disp('%')
disp('set(haxs, ''FontSize'', 14, ''Box'', ''on'', ''XGrid'', ''on'', ''YGrid'', ''on'')')

