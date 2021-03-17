classdef helperScanPatternDisplay < handle
    % This class is for plotting visualizations of a radarDataGenerator's
    % scan pattern.
    
    properties
        azLim
        elLim
        azLimTotal
        elLimTotal
        scanPtsAz
        scanPtsEl
        fov        
    end
    
    properties(Hidden)        
        % graphics handles
        animatedFig
        currentFieldIndicator
        beamFigCone
        beamFigBoresightAz
        beamFigBoresightEl
        
        % flags
        startedGif
    end
   
    methods
        
        function obj = helperScanPatternDisplay( radar )
            
            if strcmp(radar.ScanMode,'Electronic')
                obj.azLim = radar.ElectronicAzimuthLimits;
                obj.elLim = radar.ElectronicElevationLimits;
            elseif strcmp(radar.ScanMode,'Mechanical')
                obj.azLim = radar.MechanicalAzimuthLimits;
                obj.elLim = radar.MechanicalElevationLimits;
            else
               error('Combined electronic/mechanical scanning not supported'); 
            end
            
            obj.fov = radar.FieldOfView;
            
            % number of required points
            numPtsAz = floor(diff(obj.azLim)/obj.fov(1)) + 1;
            numPtsEl = floor(diff(obj.elLim)/obj.fov(2)) + 1;
            
            % amount of overscan
            osAz = numPtsAz*obj.fov(1) - diff(obj.azLim);
            osEl = numPtsEl*obj.fov(2) - diff(obj.elLim);
            
            % actual scanned area (splitting oversteer evenly)
            obj.azLimTotal = obj.azLim + [-osAz osAz]/2;
            obj.elLimTotal = obj.elLim + [-osEl osEl]/2;
            
            % contract by half-fov to get endpoints
            azScanLimTotal = obj.azLimTotal + [obj.fov(1) -obj.fov(1)]/2;
            elScanLimTotal = obj.elLimTotal + [obj.fov(2) -obj.fov(2)]/2;
            
            obj.scanPtsAz = linspace(azScanLimTotal(1), azScanLimTotal(2), numPtsAz);
            obj.scanPtsEl = linspace(elScanLimTotal(1), elScanLimTotal(2), numPtsEl);
            
        end
        
        function updatePlot( obj,ang,fname )
            
            if isempty(obj.animatedFig)
                obj.initializeAnimatedPlot;
            end
            
            subplot(1,2,1);
            
            x1 = ang(1) - obj.fov(1)/2;
            x2 = ang(1) + obj.fov(1)/2;
            y1 = ang(2) - obj.fov(2)/2;
            y2 = ang(2) + obj.fov(2)/2;
            
            obj.currentFieldIndicator.XData = [x1 x1 x2 x2];
            obj.currentFieldIndicator.YData = [y1 y2 y2 y1];
            obj.currentFieldIndicator.FaceColor = 'magenta';
            
            subplot(1,2,2);
            
            [x,y,z] = obj.drawCone( obj.fov,ang );
            obj.beamFigCone.XData = x;
            obj.beamFigCone.YData = y;
            obj.beamFigCone.ZData = z + 1;
            
            [u,v,w] = sph2cart(ang(1)*pi/180,ang(2)*pi/180,1);
            
            m = sqrt(u^2+v^2);
            obj.beamFigBoresightAz.UData = u/m;
            obj.beamFigBoresightAz.VData = v/m;
            obj.beamFigBoresightAz.WData = 0;
            
            m = sqrt(u^2+w^2);
            obj.beamFigBoresightEl.UData = u/m;
            obj.beamFigBoresightEl.VData = 0;
            obj.beamFigBoresightEl.WData = w/m;

            % rotate quiver head into XZ plane
%             vd = obj.beamFigBoresightEl.Head.VertexData;
%             vd = vd - [0;1;1];
%             vd = roty(-ang(2))*rotx(90)*roty(ang(2))*vd;
%             vd = vd + [0;1;1];
%             vd = single(vd);
%             obj.beamFigBoresightEl.Head.VertexData = vd;
            
            drawnow nocallbacks;
            
            if nargin > 2
                frame = getframe(gcf);
                im = frame2im(frame);
                [imind,cm] = rgb2ind(im,256);
                if obj.startedGif
                    imwrite(imind,cm,fname,'gif','WriteMode','append','DelayTime',0.1);
                else
                    imwrite(imind,cm,fname,'gif', 'Loopcount',inf);
                    obj.startedGif = true;
                end
            end
            
        end
        
        function makeOverviewPlot( obj )
            % Display the overview plot of the scan pattern
            
            [scanPtsAzGrid,scanPtsElGrid] = meshgrid(obj.scanPtsAz,obj.scanPtsEl);
            h(1) = plot(scanPtsAzGrid(:),scanPtsElGrid(:),'*');
            
            h(2) = line([obj.scanPtsAz(1) obj.scanPtsAz(1)]-obj.fov(1)/2,obj.elLimTotal,'color',[0 0 0 .3]);
            for ind = 1:numel(obj.scanPtsAz)
                line([obj.scanPtsAz(ind) obj.scanPtsAz(ind)]+obj.fov(1)/2,obj.elLimTotal,'color',[0 0 0 .3]);
            end
            
            line(obj.azLimTotal,[obj.scanPtsEl(1) obj.scanPtsEl(1)]-obj.fov(2)/2,'color',[0 0 0 .3]);
            for ind = 1:numel(obj.scanPtsEl)
                line(obj.azLimTotal,[obj.scanPtsEl(ind) obj.scanPtsEl(ind)]+obj.fov(2)/2,'color',[0 0 0 .3]);
            end
            
            h(3) = line([obj.azLim(1) obj.azLim(1)],[obj.elLim(1) obj.elLim(2)],'color','red');
            line([obj.azLim(2) obj.azLim(2)],[obj.elLim(1) obj.elLim(2)],'color','red');
            line([obj.azLim(1) obj.azLim(2)],[obj.elLim(1) obj.elLim(1)],'color','red');
            line([obj.azLim(1) obj.azLim(2)],[obj.elLim(2) obj.elLim(2)],'color','red');
            
            xlabel('Azimuth');
            ylabel('Elevation');
            legend(h,'Scan Points','FoV','ScanLimits','AutoUpdate','off');
            axis equal;
            xlim(obj.azLimTotal+[-obj.fov(1) obj.fov(1)]/2);
            ylim(obj.elLimTotal+[-obj.fov(2) obj.fov(2)]/2);
            title('Scan Pattern');
            
        end
        
        function initializeAnimatedPlot( obj )
            % Initialize the beam plot, do not display until update
            
            obj.animatedFig = figure('Color','white');
            obj.animatedFig.Position(3) = 2*obj.animatedFig.Position(3);
            
            subplot(1,2,1);
            obj.makeOverviewPlot;
            legend('off');
            
            % initialize current field indicator
            hold on;
            obj.currentFieldIndicator = patch('Vertices',[0 0;0 1;1 1;1 0].*obj.fov(:).','Faces',[1 2 3 4],...
                'FaceColor','none','EdgeColor','none','FaceAlpha',.2);
            hold off;
            
            subplot(1,2,2);
            
            orange = [253 106 2]/256;
            
            obj.beamFigCone = mesh(nan(2),nan(2),nan(2),'FaceColor','cyan','FaceAlpha',.5,'EdgeColor',orange,'EdgeAlpha',.2);
            hold on;
            obj.beamFigBoresightAz = quiver3(0,0,0,1,0,0,'linewidth',2,'autoscale','off','MaxHeadSize',0.75);
            obj.beamFigBoresightEl = quiver3(0,1,1,1,0,0,'linewidth',2,'autoscale','off','MaxHeadSize',0.75);
            
            line([0 0],[0 1],[1 1],'color',[0 0 0 .5],'linestyle','--');
            line([0 0],[0 0],[0 1],'color',[0 0 0 .5],'linestyle','--');
            
            % plot az/el radial grid lines
            for t = -80:10:80
                line([0 cosd(t)],[0 sind(t)],[0 0],'color',[0 0 0 0.3]);
            end
            for t = -80:10:80
                line([0 cosd(t)],[1 1],[1 1+sind(t)],'color',[0 0 0 0.3]);
            end
            
            hold off;
            xlabel('X'); ylabel('Y'); zlabel('Z');
            axis equal;
            xlim([0 1]);
            ylim([-1 1]);
            zlim([0 2]);
            
            title('Beam Pointing');
            
            set(gca,'xtick',[0 1]);
            set(gca,'ytick',[-1 0 1]);
            set(gca,'ztick',[0 1 2]);
            
        end
        
    end
    
    methods(Static)
        
        function [ X,Y,Z ] = drawCone( fov,ang )
            
            ka = tand(fov(1)/2);
            kb = tand(fov(2)/2);
            
            % Number of points along central axis and around circumference
            Nr = 10;
            Nt = 20;
            
            % base cone
            raz = linspace(0,ka,Nr);
            rel = linspace(0,kb,Nr);
            t = linspace(0,2*pi,Nt);
            Y = raz.'*cos(t);
            Z = rel.'*sin(t);
            X = sqrt(Y.^2/ka^2 + Z.^2/kb^2);
            
            R = rotz(ang(1))*roty(-ang(2));
            
            [X,Y,Z] = deal(R(1,1)*X + R(1,2)*Y + R(1,3)*Z,...
                           R(2,1)*X + R(2,2)*Y + R(2,3)*Z,...
                           R(3,1)*X + R(3,2)*Y + R(3,3)*Z);
            
        end
        
    end
    
end