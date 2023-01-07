%% Change domain color %%
    
    %% if only 1 plot
    hContour = fig.Children(2).Children(1)
    %% if subplot
    hContour = fig.Children.Children(1).Children(5);
    %% for electrode
    fig.Children(2).Children(10)
    %%
    hFills=hContour.FacePrims;
            [hFills.ColorType] = deal('truecoloralpha');  % default = 'truecolor'
            try
                hFills(1).ColorData = uint8([150;150;150;255]);
                for idx = 2 : numel(hFills)
                    hFills(idx).ColorData(4) = 0;   % default=255
                end
            catch
            end
            
    %% Potential well combined 4 plots
  for ii=2:2:8  
        hContour = fig.Children(ii).Children(1);
        hFills=hContour.FacePrims;
                [hFills.ColorType] = deal('truecoloralpha');  % default = 'truecolor'
                try
                    hFills(1).ColorData = uint8([150;150;150;255]);
                    for idx = 2 : numel(hFills)
                        hFills(idx).ColorData(4) = 0;   % default=255
                    end
                catch
                end
  end