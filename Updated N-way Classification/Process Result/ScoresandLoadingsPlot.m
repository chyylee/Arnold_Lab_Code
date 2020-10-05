%% ScoresandLoadingsPlot.m with ellipses
% Plots scores and loadings for PCA, PLSDA, and PLSR scripts (this function
% is used within PCA.m, PLSR.m, and PLSDA.m)
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Melissa Lemke, Arnold Lab, University of Michigan, Biomedical Engineering
% March 16th, 2018
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function [pep] = ScoresandLoadingsPlot(LV1_scores, LV2_scores, LV3_scores, ...
    LV1_loads, LV2_loads,LV3_loads,yblock, scores_xlabel,scores_ylabel,...
    scores_zlabel, loads_xlabel, loads_ylabel, loads_zlabel,...
    classes,  xnames, colorbarlabel, ttl, filettl, errortext, isPLSR, ...
    scores_num_LVs, loads_num_LVs)
% These workspaces contain colormaps and schemes I've made (add folder
% containing them to your path), you can develop your own, use these,
% or use those built into Matlab
% load('colormaps.mat') right now we are using Matlab's 'jet' map
load('colorschemes','quick20','pastels');
% colors = quick20([1:4 11:14 16:19],:); % this chooses the colorscheme
colors=pastels;
f1=figure();

f1.PaperUnits = 'inches';
%size for ppts, you may need to go into Matlabs 'copy options' settings to
%tell it to listen to this and not whatever the figure size is when you
%copy it
f1.PaperPosition = [0 0 9.25 4.5];
%size for papers
% f2.PaperPosition = [0 0 3 2];

% This makes a cell for each class containing all scores for the respective
% LV, not used for PLSR bc the yblock is not logical i.e. no classes
if ~isPLSR
    for i = 1:length(yblock(1,:))
        LV1{i} = LV1_scores((find(yblock(:,i)==1)));
        if scores_num_LVs>1
            LV2{i} = LV2_scores((find(yblock(:,i)==1)));
            if scores_num_LVs>2
                LV3{i} = LV3_scores((find(yblock(:,i)==1)));
            end
        end
    end
end

%% Scores plot
if scores_num_LVs < 3; subplot(10,2,[3 5 7 9 11 13 15 17]); end


if isPLSR % PLSR scores plot
    colormap(jet)
    %caxis([0.5 5]) % can fix the colorbar limits
    if scores_num_LVs == 1
        LV2_scores = 0*LV1_scores;
    elseif scores_num_LVs == 3
        scatter3(LV1_scores,LV2_scores,LV3_scores,40,yblock,'filled')
    else
        scatter(LV1_scores,LV2_scores,40,yblock,'filled')
    end
    c = colorbar;
    ylabel(c, colorbarlabel)
else % PCA and PLSDA scores plot
    for i = 1:length(yblock(1,:))
        if scores_num_LVs == 1
            LV2{i} = 0*LV1{i};
        end
        if scores_num_LVs == 3
            plot3(LV1{i},LV2{i},LV3{i},'o','MarkerFaceColor',colors(i,:),'MarkerEdgeColor','k','MarkerSize',8)
            % Create the dashed lines along the axes
            hold on
        else
            plot(LV1{i},LV2{i},'o','MarkerFaceColor',colors(i,:),'MarkerEdgeColor','k','MarkerSize',8)
            hold on
            ellip = subgroupcl([LV1{i},LV2{i}],0.9); % Adds confidence ellipses
            pep(i,:,:) = ellip.Vertices;
        end
    end
    if scores_num_LVs == 3
    xl = xlim;
            yl = ylim;
            zl = zlim;
            xlim manual
            ylim manual
            zlim manual
            if yl(1)<.01 && yl(2)>-0.1
                line(xl,[0 0],[0 0],'color','k','LineStyle','--');
            end
            if xl(1)<.01 && xl(2)>-0.1
                line([0 0],yl,[0 0],'color','k','LineStyle','--');
            end
            if zl(1)<.01 && zl(2)>-0.1
                line([0 0],[0 0],zl,'color','k','LineStyle','--');
            end
            hold on
    else
        % Create the dashed lines along the axes
            xl = xlim;
            yl = ylim;
            xlim manual
            ylim manual
            if yl(1)<.01 && yl(2)>-0.-1
                hline = line(xl,[0 0],'color','k','LineStyle','--');
            end
            if xl(1)<.01 && xl(2)>-0.-1
                vline = line([0 0],yl,'color','k','LineStyle','--');
            end
            hold on
    end
    legend(classes,'Location','northeast','Color','none','AutoUpdate','off')
end



% Title and axes labels
title('Scores Plot')
box on
xlabel(scores_xlabel)
if scores_num_LVs > 1
    ylabel(scores_ylabel)
end
if scores_num_LVs > 2
    zlabel(scores_zlabel) 
end
set(gca,'TickLabelInterpreter','none','FontSize',8)
hold off

%% Loadings plot
if loads_num_LVs < 3; subplot(10,2,[4 6 8 10 12 14 16 18]); end
if loads_num_LVs ==3; figure(); end
colors = {[.3 .3 .3]};

if loads_num_LVs==1 % bar graph if only 1D
    [loadings, index] = sort(LV1_loads); % sorted so they look nice
    barh(loadings,'FaceColor',colors{1})
    set(gca, 'ytick',1:length(loadings))
    ylim([0.5 length(loadings)+0.5])
    yticklabels(xnames(index));
    title('Loadings Plot')
    xlabel(loads_xlabel)
    set(gca,'TickLabelInterpreter','none','FontSize',8)
    propertyeditor('on') % this toggle allows the figure to open in a state
    propertyeditor('off') % where you can move things like the legend
    
elseif loads_num_LVs == 2 % Scatter plot
    hold on
    box on
    scatter(LV1_loads,LV2_loads,8,'filled','markerfacecolor',colors{1},...
        'markeredgecolor','k')
    text(LV1_loads,LV2_loads, xnames','FontSize',6,'Interpreter','none') % factor labels
    
    % Create the dashed lines along the axes
    xl = xlim;
    yl = ylim;
    xlim manual
    ylim manual
    if yl(1)<.01 && yl(2)>-0.1
        line(xl,[0 0],'color','k','LineStyle','--');
    end
    if xl(1)<.01 && xl(2)>-0.1
        line([0 0],yl,'color','k','LineStyle','--');
    end
    xlabel(loads_xlabel)
ylabel(loads_ylabel)
zlabel(loads_zlabel)
elseif loads_num_LVs == 3
    plot3(LV1_loads,LV2_loads,LV3_loads,'o','MarkerFaceColor',...
        colors{1,:},'MarkerEdgeColor','k','MarkerSize',8)
    text(LV1_loads,LV2_loads,LV3_loads, xnames','FontSize',6,'Interpreter','none') % factor labels
        box on
    % Create the dashed lines along the axes
    xl = xlim;
    yl = ylim;
    zl = zlim;
    xlim manual
    ylim manual
    zlim manual
    if yl(1)<.01 && yl(2)>-0.1
        line(xl,[0 0],[0 0],'color','k','LineStyle','--');
    end
    if xl(1)<.01 && xl(2)>-0.1
        line([0 0],yl,[0 0],'color','k','LineStyle','--');
    end
    if zl(1)<.01 && zl(2)>-0.1
        line([0 0],[0 0],zl,'color','k','LineStyle','--');
    end
    xlabel(loads_xlabel)
ylabel(loads_ylabel)
zlabel(loads_zlabel)
end

title('Loadings Plot')

set(gca,'TickLabelInterpreter','none','FontSize',8)
% propertyeditor('on') % this toggle allows the figure to open in a state
% propertyeditor('off') % where you can move things like the legend
hold off


%fake axis to get a stable position for errors and big title
a = axes;
%// Set the title and get the handle to it
ht = title(ttl,'FontSize',14);
%// Turn the visibility of the axes off
a.Visible = 'off';
%// Turn the visibility of the title on
ht.Visible = 'on';

% text for the botton of the figure, typically the cal and cv errors in
% PLSDA
if ~isempty(errortext)
    xPos = 0.48;
    yPos = -0.07;
    t=text(xPos,yPos,errortext,'Parent',a,'FontSize',8);
    set(t, 'HorizontalAlignment', 'center');
end

%saves the figure
savefig(f1,strcat(datestr(today()),strtok(extractAfter(filettl,11),'.'),'.fig'))

end