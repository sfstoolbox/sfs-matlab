function [tau, R, Delta] = retarded_time(x, t, vs, conf)
    % x  [nx3]
    % t  [nx1]
    % vs [1x3] or [nx3]

    c = conf.c;
    
    v = vector_norm(vs,2);  % velocity of sound source
    nvs = bsxfun(@rdivide,vs,v);  % direction of movement
    M = v/c;  % Mach number
    
    r = vector_norm(bsxfun(@minus,x,bsxfun(@times,t,vs)),2);  % r = |x0-xs| 
    % component of x - vs*t is the direction of nvs
    xparallel = sum( bsxfun(@times, x, nvs), 2) - t.*v;  
    Delta = sqrt( M.^2.*xparallel.^2 + (1-M.^2).*(r.^2) );
    R = (M.*xparallel + Delta)./(1-M.^2);
    tau = R./c;
end