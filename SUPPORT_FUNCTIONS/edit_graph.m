% Function developed by Mike Garrity Dec. 31, 2015
% https://blogs.mathworks.com/graphics/2015/12/31/interactive-graph-layout/


function edit_graph(f,h)
    
    % Figure out which node is closest to the mouse. 
    a = ancestor(h,'axes');
    pt = a.CurrentPoint(1,1:2);
    dx = h.XData - pt(1);
    dy = h.YData - pt(2);
    len = sqrt(dx.^2 + dy.^2);
    [lmin,idx] = min(len);
    
    % If we're too far from a node, just return
    tol = max(diff(a.XLim),diff(a.YLim))/20;    
    if lmin > tol || isempty(idx)
        return
    end
    node = idx(1);

    % Install new callbacks on figure
    f.WindowButtonMotionFcn = @motion_fcn;
    f.WindowButtonUpFcn = @release_fcn;

    % A ButtonMotionFcn that changes XData & YData
    function motion_fcn(~,~)
        newx = a.CurrentPoint(1,1);
        newy = a.CurrentPoint(1,2);
        h.XData(node) = newx;
        h.YData(node) = newy;
        drawnow;
    end

    % A ButtonUpFcn which stops dragging
    function release_fcn(~,~)
        f.WindowButtonMotionFcn = [];
        f.WindowButtonUpFcn = [];
    end
end

