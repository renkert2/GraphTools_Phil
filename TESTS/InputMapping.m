cmap = [
    2 3;
    3 1;
    5 4;
    4 6
    6 2];

for i = 1:6
    inputs(i) = GraphInput(sprintf("Input %d",i));
end

input_conn_map = arrayfun(@(x) inputs(x), cmap);


cmap = ReformatSerialConnections(cmap)
input_conn_map = ReformatSerialConnections(input_conn_map)

function conn_map = ReformatSerialConnections(conn_map)
inters = intersect(conn_map(:,1),conn_map(:,2)); % Find common elements in first and second columns of conn_map
if ~isempty(inters) % If intersections exist
    for i = 1:numel(inters)
        target_val = conn_map(inters(i) == conn_map(:,1),2);
        conn_map(inters(i) == conn_map(:,2),2) = target_val;
    end
end
end