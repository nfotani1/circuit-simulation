function bool = creates_cycle (reset, A_column)
% Call this routine as each edge added to a graph.
% Whenever an edge creates a cycle with previously added
% edges, this routine returns "true"; otherwise, it
% returns "false".
% Call the routine with reset = true to start over
% by removing all edges.

persistent parent;

if (nargin==1 && reset); bool = 0; parent = []; return; end

if (reset || isempty(parent))
    parent = - ones(length(A_column),1);
end

vertex_of_edge = find(A_column~=0);

if (length(vertex_of_edge)~=2)
    disp('A_column is flawed. It doesn''t have exactly 2 nonzero entries');
    return;
end
src = vertex_of_edge(1); dst = vertex_of_edge(2);

parent_of_src = find_parent(src);
parent_of_dst = find_parent(dst);

 if (parent_of_src==parent_of_dst)
    bool = true;
else
    bool = false;
    parent_of_parent_of_src = find_parent(parent_of_src);
    parent_of_parent_of_dst = find_parent(parent_of_dst);
    parent(parent_of_parent_of_src) = parent_of_parent_of_dst;
end

% fprintf('%2i ',1:length(A_column)); fprintf('\n');
% fprintf('%2i ',parent); fprintf('\n');

    function parent_of_p = find_parent(p)
        p_temp = parent(p);
        while (p_temp ~= -1)
            p = p_temp;
            p_temp = parent(p);
        end
        parent_of_p = p;
    end

end