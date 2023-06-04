function h=display_model(model)
    hold on;
    h=plot_polygon(model.walls);
    for i=1:length(h)
        set(h(i),'color','red');
        set(h(i),'linewidth',2);
    end

    for i=1:length(model.reflectors)
        plot(model.reflectors{i}(:,1),model.reflectors{i}(:,2),'g','linewidth',2);
    end

    for i=1:length(model.obstacles)
        h=plot_polygon(model.obstacles{i});
        for j=1:length(h)
            set(h(j),'color','blue');
            set(h(j),'linewidth',2);
        end
    end
end