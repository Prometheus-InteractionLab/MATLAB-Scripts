function ray = createray(startpos,mode,vararg,power,indexOfRefraction)
% Create a ray starting at startpos and with a direction defined either by
% the user or by a second point along the ray
%startpos: [x0;y0]
%mode: 'dir' or 'point'
%vararg: [u,v] or [x1;y1]
%power: ray power
%indexOfRefraction: index of refraction for the material that the ray is
%passing through

%Define unit vector for ray
switch(mode)
    case 'dir'
        %user inputs a direction for vararg
        normdir=normc(vararg);
    case 'point'
        %user inputs the second point for vararg
        normdir=normc(vararg-startpos);
    otherwise
        error('mode does not exist')
        beep
        return
end

%ray has a starting point
ray.start=startpos;
%parallel unit vector for the ray
ray.pu=normdir;
%function reference - gives position as a function of parameter t
t_ray=@ (t) startpos + t*normdir;
ray.t_ray=t_ray;
ray.power=power;
ray.n=indexOfRefraction;

end