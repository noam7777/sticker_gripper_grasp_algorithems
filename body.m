classdef body %required to send reference to the class methods and to use recursive functions
    %BODY - a class of a closed 2D body
    %   it has n points and n segments because it is closed
    
    properties
        n % the number of points on the body, without cycling over itselves
        Cx % center of mass 
        Cy
        per % the entire perimiter of the body
        psi     % [x,y,s] coordinates of the body,s is the total track from point 1 to the current point s(i)=s(i-1)+ds(i-1)
        dpsi    % differance array of psi, contains 1 less element, psi(:,3) is the track from the current point to the next    
        Tan     % an array of unit vectors that represent the tangent directions of the body segments
    end
    
    methods
        function obj = body(x,y,Cx,Cy)
            %BODY Construct an instance of this class
            %   Detailed explanation goes here
            obj.Cx = Cx;    % defining center of mass
            obj.Cy = Cy;
            obj.n = length(x);    % the amount of body points (and colsed body segments)      
            n = obj.n;
            obj.psi = zeros(n*3+1,3);           % in case of a closed body, a periodic definition of the body
                                                % points is required, 3 cycles is enough. the two edges of the
                                                % point array repesents the same body point hence the +1 in
                                                % the psi array
            obj.dpsi = zeros(n*3,3);    % diffrence array, point(i+1) - point(i), the size is the 
                                        % amount of body segments
            obj.Tan = zeros(n*3,2);

            
            obj.psi(:,1) = [x,x,x,x(1)];
            obj.psi(:,2) = [y,y,y,y(1)];
            obj.dpsi(:,1:2) = obj.psi(2:end,1:2) - obj.psi(1:end-1,1:2);
            obj.dpsi(:,3) = sqrt(sum(obj.dpsi(:,1:2).^2,2));
            
            obj.psi(1,3)=0;
            for i = 2:3*obj.n+1
                obj.psi(i,3) = obj.psi(i-1,3)+obj.dpsi(i-1,3);
            end
            obj.per = obj.psi(n+1,3); 
            obj.Tan = obj.dpsi(:,1:2)./obj.dpsi(:,3);
            
        end
        function print(obj) 
            plot(obj.psi(1:obj.n+1,1),obj.psi(1:obj.n+1,2))
            hold on
            scatter(obj.Cx,obj.Cy,100,'k')
            scatter(obj.Cx,obj.Cy,100,'+k')
            hold off
        end
        function test_body_shape(obj,dt) 
            plot(obj.psi(1:obj.n+1,1),obj.psi(1:obj.n+1,2))
            hold on
            scatter(obj.Cx,obj.Cy,100,'k')
            scatter(obj.Cx,obj.Cy,100,'+k')
            for ii = 1:obj.n
                scatter(obj.psi(ii,1),obj.psi(ii,2),5)
                pause(dt)
            end
            hold off
        end
        function plot_body_point(obj,sO)
            xO = interp1(obj.psi(:,3),obj.psi(:,1),sO);   % x coordinates of E
            yO = interp1(obj.psi(:,3),obj.psi(:,2),sO);   % y coordinates of E
                scatter(xO,yO,30,'k')
            xO2 = interp1(obj.psi(:,3),obj.psi(:,1),sO+obj.per);   % x coordinates of E
            yO2 = interp1(obj.psi(:,3),obj.psi(:,2),sO+obj.per);   % y coordinates of E
                scatter(xO2,yO2,20,'k')
        end
    end
end

