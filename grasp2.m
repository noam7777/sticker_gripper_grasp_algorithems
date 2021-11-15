classdef grasp2
    %GRASP2 a representation of a two film grasp
    %   Detailed explanation goes here
    
    properties
        body
        filmA
        filmB
        sO
        sA
        sB
        SA
        SB
        LA
        LB
        A
        B
        d
        LA_loose_cond_geo
        LB_loose_cond_geo
        Ax
        Ay
        Bx
        By
        Dx
        Dy
        
        ADx_unity
        ADy_unity
        BDx_unity
        BDy_unity
        CDx_unity
        CDy_unity
        CBx_unity
        CBy_unity
        CAx_unity
        CAy_unity
        
        TA_unity
        TB_unity
        
        Dx_mf
        Dy_mf
        LA_loose_cond
        LB_loose_cond
        
        ADx_unity_mf
        ADy_unity_mf
        BDx_unity_mf
        BDy_unity_mf
        CDx_unity_mf
        CDy_unity_mf
        
        TanA_x
        TanA_y
        TanB_x
        TanB_y
        
        alfa
        beta
        
        theta1
        theta2
        
        TAcr
        TBcr
        
        mg_cr_A
        mg_cr_B
        mg_cr
    end
    
    methods
        function obj = grasp2(body,sO,filmA,filmB,s_tot)
            %GRASP2 Construct an instance of this class
            %   Detailed explanation goes here
            obj.body = body;
            obj.sO = sO;
            obj.filmA = filmA;  %the relative progerss of the peeling interface
            obj.filmB = filmB;  
            nA = filmA.n;
            nB = filmB.n;
            %% determinning the sample points in the sA sB plane
            sA = linspace(0,s_tot,filmA.n);
            sB = linspace(0,s_tot,filmB.n);
            [SB,SA] = meshgrid(sB,sA);
            obj.sA = sA;
            obj.sB = sB;
            obj.SA = SA;
            obj.SB = SB;           sA_actual = sO + body.per + sA; % the actual position of the film on the body, between the two ofset zones
            sB_actual = sO + body.per - sB;
            [SB_actual,SA_actual] = meshgrid(sB_actual,sA_actual);
            %% calculating the free film length 
            lA = sA+filmA.l_init;
            lB = sB+filmB.l_init;
            [LB,LA] = meshgrid(lB,lA);   
            obj.LA = LA;
            obj.LB = LB;
            %% calculating the pose of the adhesion interface
            A = zeros(filmA.n,2);   % array of all the possible A positions
            B = zeros(filmB.n,2);
            A(:,1) = interp1(body.psi(:,3),body.psi(:,1),sA_actual);   % x coordinates of A
            A(:,2) = interp1(body.psi(:,3),body.psi(:,2),sA_actual);   % y coordinates of A
            B(:,1) = interp1(body.psi(:,3),body.psi(:,1),sB_actual);   % x coordinates of A
            B(:,2) = interp1(body.psi(:,3),body.psi(:,2),sB_actual);   % y coordinates of A
            obj.A = A;
            obj.B = B;
            [Bx,Ax] = meshgrid(obj.B(:,1),obj.A(:,1));
            [By,Ay] = meshgrid(obj.B(:,2),obj.A(:,2)); %turning the points to an (sA,sB) matrices - sA is the row , sB is colomn
            obj.Ax = Ax;
            obj.Ay = Ay;
            obj.Bx = Bx;
            obj.By = By;
            d = sqrt((Ax-Bx).^2+(Ay-By).^2);    % the distance between A and B
            obj.d=d;
        
            
            
            %% loose films check
            % the s_index positions in which one of the filmes is loose
            % for geometrical resons, where one circle contain the other.
            LA_loose_cond_geo = LA>=(LB+d-10*eps); % A is loose. consider removing equal sign.
            LB_loose_cond_geo = LB>=(LA+d-10*eps); % B is loose
            L_loose_cond_geo = LA_loose_cond_geo | LB_loose_cond_geo; % this condition do not include the sA = sB = 0 case
            obj.LA_loose_cond_geo = LA_loose_cond_geo;
            obj.LB_loose_cond_geo = LB_loose_cond_geo;
            if max(max(LA+LB<d)) % this infers a serious geometrical error
                disp("geometrical error! the sum of lA and lB is smaller then the distance between the A and B points, the circles are strangers")
                return;
            end
            
            %% calculating the position of the griper, position D
%             Dx = zeros(nA,nB);
%             Dy = zeros(nA,nB);
%             [sA_ind,sB_ind] = find(~L_loose_cond_geo); %temp
            
            a = (d.^2+LA.^2-LB.^2)./(2.*d);
            h = (LA.^2-a.^2).^0.5;
            p2x = Ax +a./d.*(Bx-Ax);
            p2y = Ay +a./d.*(By-Ay);
            Dx = zeros(nA,nB);
            Dy = zeros(nA,nB);
            %calculating all the simple elements of D, where the tow
            %circles intercect definetly. the other elements remain zero.
            Dx(~L_loose_cond_geo) = p2x(~L_loose_cond_geo) - h(~L_loose_cond_geo).*(By(~L_loose_cond_geo)-Ay(~L_loose_cond_geo))./d(~L_loose_cond_geo);
            Dy(~L_loose_cond_geo) = p2y(~L_loose_cond_geo) + h(~L_loose_cond_geo).*(Bx(~L_loose_cond_geo)-Ax(~L_loose_cond_geo))./d(~L_loose_cond_geo);
            obj.Dx = Dx;
            obj.Dy = Dy;
            
            %spacial case1 (lA < lB+d) && (lB < lA + d) will be calculated
            %later, after the mass distribution factor will be added to the
            %loose conditon
            
            %% calculating the unit vectors from the interface points and the center of mass to the position of the griper
            ADx = obj.Dx-Ax; 
            ADy = obj.Dy-Ay;
            AD_abs = (ADx.^2+ADy.^2).^0.5;
            ADx_unity = ADx./AD_abs;
            ADy_unity = ADy./AD_abs;    
            obj.ADx_unity = ADx_unity; % checked
            obj.ADy_unity = ADy_unity;

            BDx = obj.Dx-Bx;
            BDy = obj.Dy-By;
            BD_abs = (BDx.^2+BDy.^2).^0.5;
            BDx_unity = BDx./BD_abs;
            BDy_unity = BDy./BD_abs; 
            obj.BDx_unity = BDx_unity; % checked
            obj.BDy_unity = BDy_unity;
            
            CDx = obj.Dx-body.Cx;
            CDy = obj.Dy-body.Cy;
            CD_abs = (CDx.^2+CDy.^2).^0.5;
            CDx_unity = CDx./CD_abs;
            CDy_unity = CDy./CD_abs;
            obj.CDx_unity = CDx_unity; % checked
            obj.CDy_unity = CDy_unity;        
            
            CBx = Bx - obj.body.Cx;
            CBy = By - obj.body.Cy;
            CB_abs = (CBx.^2+CBy.^2).^0.5;
            CBx_unity = CBx./CB_abs;
            CBy_unity = CBy./CB_abs;            
            obj.CBx_unity = CBx_unity;
            obj.CBy_unity = CBy_unity;
            
            CAx = Ax - obj.body.Cx;
            CAy = Ay - obj.body.Cy;
            CA_abs = (CAx.^2+CAy.^2).^0.5;
            CAx_unity = CAx./CA_abs;
            CAy_unity = CAy./CA_abs;
            obj.CAx_unity = CAx_unity;
            obj.CAy_unity = CAy_unity;
            %% definining D at start position
            if sA(1) == 0 && sB(1) == 0
                obj.Dx(1,1) = CAx_unity(1).*filmA.l_init+Ax(1);
                obj.Dy(1,1) = CAy_unity(1).*filmA.l_init+Ay(1);
            end            
        
            %%  calculate the tention destribution between the two films
            
            det = ADx_unity.*BDy_unity-BDx_unity.*ADy_unity; % for cramer solution for linear system
            
            
            TA_unity = (CDx_unity.*BDy_unity-BDx_unity.*CDy_unity)./det;
            TB_unity = (CDy_unity.*ADx_unity-ADy_unity.*CDx_unity)./det;

            
            %% if TA or TB is smaller then 0, it is considered loose and D,T,(alfa i not yet setteled) need to be refreshed
            %% recalculating D and TA,TB when one of the films is loosed
            Dx_mf = Dx;         
            Dy_mf = Dy;
            %for a loosed film_A:
            LA_loose_cond_mass = TA_unity<=0;
            LA_loose_cond = LA_loose_cond_mass|LA_loose_cond_geo;
            Dx_mf(LA_loose_cond) = Bx(LA_loose_cond)+CBx_unity(LA_loose_cond).*LB(LA_loose_cond);
            Dy_mf(LA_loose_cond) = By(LA_loose_cond)+CBy_unity(LA_loose_cond).*LB(LA_loose_cond);
            TA_unity(LA_loose_cond) = 0;
            TB_unity(LA_loose_cond) = 1;

            %for a loosed film_B:
            LB_loose_cond_mass = TB_unity<=0;
            LB_loose_cond = LB_loose_cond_mass|LB_loose_cond_geo;
            Dx_mf(LB_loose_cond) = Ax(LB_loose_cond)+CAx_unity(LB_loose_cond).*LA(LB_loose_cond);
            Dy_mf(LB_loose_cond) = Ay(LB_loose_cond)+CAy_unity(LB_loose_cond).*LA(LB_loose_cond);
            TA_unity(LB_loose_cond) = 1;
            TB_unity(LB_loose_cond) = 0;   
            obj.TA_unity = TA_unity;
            obj.TB_unity = TB_unity;
            
            obj.Dx_mf = Dx_mf;
            obj.Dy_mf = Dy_mf;
            obj.LA_loose_cond = LA_loose_cond;
            obj.LB_loose_cond = LB_loose_cond;
            %% now we fix AD,BD where needed
            ADx_mf = Dx_mf-Ax; 
            ADy_mf = Dy_mf-Ay;
            AD_abs_mf = (ADx_mf.^2+ADy_mf.^2).^0.5;
            ADx_unity_mf = ADx_mf./AD_abs_mf;
            ADy_unity_mf = ADy_mf./AD_abs_mf;    
            obj.ADx_unity_mf = ADx_unity_mf; % checked
            obj.ADy_unity_mf = ADy_unity_mf;

            BDx_mf = Dx_mf-Bx;
            BDy_mf = Dy_mf-By;
            BD_abs_mf = (BDx_mf.^2+BDy_mf.^2).^0.5;
            BDx_unity_mf = BDx_mf./BD_abs_mf;
            BDy_unity_mf = BDy_mf./BD_abs_mf; 
            obj.BDx_unity_mf = BDx_unity_mf; % checked
            obj.BDy_unity_mf = BDy_unity_mf;
            
            CDx_mf = obj.Dx_mf-body.Cx;
            CDy_mf = obj.Dy_mf-body.Cy;
            CD_abs_mf = (CDx_mf.^2+CDy_mf.^2).^0.5;
            CDx_unity_mf = CDx_mf./CD_abs_mf;
            CDy_unity_mf = CDy_mf./CD_abs_mf;
            obj.CDx_unity_mf = CDx_unity_mf; % checked
            obj.CDy_unity_mf = CDy_unity_mf;                  
            
            
            
            
            %%  calc tangent unit vector at descrete adhesion points
            TanA_x = interp1(body.psi(1:end-1,3),body.Tan(:,1),SA_actual);
            TanA_y = interp1(body.psi(1:end-1,3),body.Tan(:,2),SA_actual);
            TanB_x = interp1(body.psi(1:end-1,3),body.Tan(:,1),SB_actual);
            TanB_y = interp1(body.psi(1:end-1,3),body.Tan(:,2),SB_actual);
            obj.TanA_x=TanA_x;
            obj.TanA_y=TanA_y;
            obj.TanB_x=TanB_x;
            obj.TanB_y=TanB_y;
            
            %% calc alfa and beta
            %alfa = zeros(filmA.n,filmB.n);
            %beta = zeros(filmA.n,filmB.n);
            x1 = -TanA_x;
            y1 = -TanA_y;
            x2 = ADx_unity_mf;
            y2 = ADy_unity_mf;
            alfa = atan2(x1.*y2-y1.*x2,x1.*x2+y1.*y2);
            
            x1 = BDx_unity_mf;
            y1 = BDy_unity_mf;            
            x2 = TanB_x;
            y2 = TanB_y;
            beta = atan2(x1.*y2-y1.*x2,x1.*x2+y1.*y2);

            % defining alfa, beta of loose films is loose
            alfa(LA_loose_cond) = nan; %the filme is loose, there is no sense at difining its angle!
            beta(LB_loose_cond) = nan;
            obj.alfa = alfa;
            obj.beta = beta;
            
            theta1 = (pi/2-alfa);
            theta2 = (pi/2-beta);
            theta1 = max(theta1,0); %TODO - what is the peeling force for "small" by the 180 angle paper? for now it is approximated to the 90deg peeling force
            theta2 = max(theta2,0);
            
            obj.theta1 = theta1;
            obj.theta2 = theta2;
            
            TAcr=1./(1-cos(alfa));
            obj.TAcr = TAcr;
            TBcr=1./(1-cos(beta));
            obj.TBcr = TBcr;
            
            mg_cr_A = TAcr./TA_unity;
            mg_cr_B = TBcr./TB_unity;
            mgcr = min(mg_cr_A,mg_cr_B);
            obj.mg_cr_A = mg_cr_A;
            obj.mg_cr_B = mg_cr_B;
            obj.mg_cr = mgcr;
            %calculate the peeling direction for each point on the plane
%             [mg_cr_gr2,pealing_index_max_gr2] = max(interp2(gr2.SB,gr2.SA,gr2.mg_cr,strm(:,1),strm(:,2)));
%             [mg_cr_gr2,pealing_index_max_gr2] = max(interp2(mesh of x axis,mesh of y axis,scalar field, point x of desired scalar,point y of desierd scalar));
            
            
        end
        function plot_sticker_at_state(obj,sa_k,sb_k)
            Ax_k = interp2(obj.SB,obj.SA,obj.Ax,sb_k,sa_k);   % x coordinates of A at the current sa,sb state
            Ay_k = interp2(obj.SB,obj.SA,obj.Ay,sb_k,sa_k);   
            Bx_k = interp2(obj.SB,obj.SA,obj.Bx,sb_k,sa_k);   
            By_k = interp2(obj.SB,obj.SA,obj.By,sb_k,sa_k);   
            Dx_k = interp2(obj.SB,obj.SA,obj.Dx,sb_k,sa_k);   
            Dy_k = interp2(obj.SB,obj.SA,obj.Dy,sb_k,sa_k);  
            plot([Ax_k,Dx_k,Bx_k],[Ay_k,Dy_k,By_k],'k')
        end
        function plot_sticker_at_state_thick(obj,sa_k,sb_k)
            Ax_k = interp2(obj.SB,obj.SA,obj.Ax,sb_k,sa_k);   % x coordinates of A at the current sa,sb state
            Ay_k = interp2(obj.SB,obj.SA,obj.Ay,sb_k,sa_k);   
            Bx_k = interp2(obj.SB,obj.SA,obj.Bx,sb_k,sa_k);   
            By_k = interp2(obj.SB,obj.SA,obj.By,sb_k,sa_k);   
            Dx_k = interp2(obj.SB,obj.SA,obj.Dx,sb_k,sa_k);   
            Dy_k = interp2(obj.SB,obj.SA,obj.Dy,sb_k,sa_k);  
            plot([Ax_k,Dx_k,Bx_k],[Ay_k,Dy_k,By_k],'m-*','linewidth',3,'MarkerSize',10)
        end        
        
        
    end
end

