classdef grasp_single_film
    %GRASP_SINGLE_FILM Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        body
        film
        sO
        sE
        L
        E
        Ex
        Ey
        Dx
        Dy
        
        EDx_unity
        EDy_unity
        
        TE_unity
        
        TanE_x
        TanE_y
        
        eta
        
        TEcr
        
        mg_cr
        
    end
    
    methods
        function obj = grasp_single_film(body,sO,filmE,film_direction,s_tot)
            %GRASP_SINGLE_FILM Construct an instance of this class
            % film_direction=1 if it peels to the positive direction (usually counter clockwise) 
            % of the curve. 
            % direction=1 for point A 
            % direction=-1 for point B

            obj.body = body;
            obj.sO = sO;
            obj.film = filmE;  %the relative progerss of the peeling interface
            %nE = filmE.n;
            %% determinning the sample points in the sE vector
            sE = linspace(0,s_tot,filmE.n);
            obj.sE = sE;
            sE_actual = sO + body.per + sE.*film_direction; % the actual position of the film on the body, between the two ofset zones
            %% calculating the free film length 
            LE = sE+filmE.l_init;
            %% calculating the pose of the adhesion interface
            E = zeros(filmE.n,2);   % array of all the possible E positions
            E(:,1) = interp1(body.psi(:,3),body.psi(:,1),sE_actual);   % x coordinates of E
            E(:,2) = interp1(body.psi(:,3),body.psi(:,2),sE_actual);   % y coordinates of E
            obj.E = E;
            Ex = E(:,1);
            Ey = E(:,2);
            obj.Ex = Ex;
            obj.Ey = Ey; 
            
            %% calculate CE and D, simular to a loose case of a double peeling case
            CEx = Ex - obj.body.Cx;
            CEy = Ey - obj.body.Cy;
            CE_abs = (CEx.^2+CEy.^2).^0.5;
            CEx_unity = CEx./CE_abs;
            CEy_unity = CEy./CE_abs;
            obj.EDx_unity = CEx_unity;
            obj.EDx_unity = CEy_unity;
            
            Dx = Ex+CEx_unity.*LE'; %CE_unity == CD_unity for the single filme case
            Dy = Ey+CEy_unity.*LE';
            TE_unity = 1;
            
            obj.Dx = Dx;
            obj.Dy = Dy;

            
            %%  calc tangent unit vector at descrete adhesion points
            TanE_x = interp1(body.psi(1:end-1,3),body.Tan(:,1),sE_actual)';
            TanE_y = interp1(body.psi(1:end-1,3),body.Tan(:,2),sE_actual)';

            obj.TanE_x=TanE_x;
            obj.TanE_y=TanE_y;
            
%% calc eta
            %alfa = zeros(filmA.n,filmB.n);
            if film_direction == 1
                x1 = -TanE_x;
                y1 = -TanE_y;
                x2 = CEx_unity;
                y2 = CEy_unity;
                eta = atan2(x1.*y2-y1.*x2,x1.*x2+y1.*y2);
            else
                x1 = CEx_unity;
                y1 = CEy_unity;            
                x2 = TanE_x;
                y2 = TanE_y;
                eta = atan2(x1.*y2-y1.*x2,x1.*x2+y1.*y2);
            end
    
            %theta = (pi/2-eta);
            %obj.theta = theta;
            obj.eta = eta;
            TEcr=1./(1-cos(eta));
            obj.TEcr = TEcr;

            mg_cr = TEcr./TE_unity;
            obj.mg_cr = mg_cr;


            
            

        end
        
        function plot_sticker_at_state(obj,sE_k)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            Ex_k = interp1(obj.sE,obj.Ex,sE_k);   % x coordinates of E
            Ey_k = interp1(obj.sE,obj.Ey,sE_k);   % y coordinates of E
            Dx_k = interp1(obj.sE,obj.Dx,sE_k);   % x coordinates of D
            Dy_k = interp1(obj.sE,obj.Dy,sE_k);   % y coordinates of D
            plot([Ex_k,Dx_k],[Ey_k,Dy_k],'Color',[0.5,0.5,0.5])
        end
        function plot_sticker_at_state_thick(obj,sE_k)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            Ex_k = interp1(obj.sE,obj.Ex,sE_k);   % x coordinates of E
            Ey_k = interp1(obj.sE,obj.Ey,sE_k);   % y coordinates of E
            Dx_k = interp1(obj.sE,obj.Dx,sE_k);   % x coordinates of D
            Dy_k = interp1(obj.sE,obj.Dy,sE_k);   % y coordinates of D
            plot([Ex_k,Dx_k],[Ey_k,Dy_k],'m-*','linewidth',3,'MarkerSize',10)
        end
    end
end

