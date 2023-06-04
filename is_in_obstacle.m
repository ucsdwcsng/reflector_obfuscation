function k = is_in_obstacle(obstacles, loc)

n_obstacles = length(obstacles);
k = false;
for i = 1:n_obstacles
    if(loc(1)>min(obstacles{i}(:,1)) && loc(1)<max(obstacles{i}(:,1)) && loc(2)>min(obstacles{i}(:,2)) && loc(2)<max(obstacles{i}(:,2)))
        k = true;
    end
    
end