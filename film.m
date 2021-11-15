classdef film < handle
    %FILM Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        s       % array of all descrete piled lengths of the film
        n       % the size of s, repesent the amunt of samples of s over the continuse peeling process
        s_tot   % the initial length of the film-body interface
        l       % free wing length
        l_init  % the initial length of the free wing
        T       % tention
        gama 	% adhesion energy
        track       % the [x y s] coordinates that define the film. the 1st point
                    % repesents the D point that is the sticker-griper griping point.
                    % so s(0) is nan because it is not on the body.
        
    end
    
    methods
        function obj = film(l_init,s_tot,gama,n_film)
            obj.l_init = l_init; 
            obj.s_tot = s_tot;
            obj.gama = gama;
            obj.n = n_film;
            %obj.l = l_init;
        end        
    end
end

