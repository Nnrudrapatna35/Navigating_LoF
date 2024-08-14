function dist = point_to_line(pt, e1, e2)
    dist_e1 = sqrt(sum((pt-e1).^2));
    dist_e2 = sqrt(sum((pt-e2).^2));
    l = sqrt(sum((e1-e2).^2));
    
    a = e1 - e2;
    b = pt - e2;
    dist_l = norm(cross(a,b)) / norm(a);
    
    if (max(dist_e1^2, dist_e2^2) - dist_l^2 > l^2)
        dist = min(dist_e1, dist_e2);
    else
        dist = dist_l;
    end
end
