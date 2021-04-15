classdef GeoVector < util.math.Vector3
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        referenceFrame = 'ITRF2005';
        referenceEllipsoid = referenceEllipsoid('Geodetic Reference System 1980');
        epoch = util.lang.Time(datetime(2005,1,1));
    end
    
    methods
        function obj = GeoVector(varargin)
            if nargin>=3
                super_args = varargin(1:3);
            end
            if nargin == 1
                super_args = varargin(1);
            end
            if nargin == 0
                super_args = {};
            end
            obj@util.math.Vector3(super_args{:});
            obj.referenceFrame = 'ITRF2005';
            if nargin == 4
                obj.referenceFrame = varargin{4};
            end
        end
    end
    
    methods       
        function ret = getEllipsoidNormal(obj)
            [lat,lon,~] = ecef2geodetic(obj.referenceEllipsoid,obj.x,obj.y,obj.z);
            [x,y,z] = enu2ecefv(0,0,1,lat,lon);
            ret = util.math.Vector3(x,y,z);
        end
        
        function ret = getPointing(obj,azimuth,elevation)
            [lat,lon,~] = ecef2geodetic(obj.referenceEllipsoid,obj.x,obj.y,obj.z);
            u = sind(elevation);
            coselev = cosd(elevation);
            n = coselev .* cosd(azimuth);
            e = coselev .* sind(azimuth);
            [x,y,z] = enu2ecefv(e,n,u,lat,lon);
            ret = util.math.Vector3(x,y,z);
        end
        
        function ret = getNedToEcef(obj)
            [lat,lon,~] = ecef2geodetic(obj.referenceEllipsoid,obj.x,obj.y,obj.z);
            ret = util.math.Quaternion.fromNed2Ecef(lat,lon);
        end
        
        function obj = lookAtSpheroid(obj,lookDirection,height)
            if nargin<3
                height=0;
            end
            a = obj.referenceEllipsoid.SemimajorAxis+height;
            b = obj.referenceEllipsoid.SemiminorAxis+height;
            x0=obj.x./a;
            y0=obj.y./a;
            z0=obj.z./b;
            ux=lookDirection.x;
            uy=lookDirection.y;
            uz=(a./b).*lookDirection.z;
            
            A = ux.^2 + uy.^2 + uz.^2;
            p2 = (ux.*x0 + uy.*y0 + uz.*z0)./A;
            q = (x0.^2+y0.^2+z0.^2-1)./A;
            t = -p2-sqrt(p2.^2-q);
            
            obj= util.math.GeoVector(a.*(x0+ux.*t),a.*(y0+uy.*t),b.*(z0+uz.*t));
            
        end
        
        function ll=plotKml(obj,varargin)
           ll = kml([char(java.util.UUID.randomUUID.toString),'.kml']);
           [lat,lon,height] = toWgs84Egm96(obj);
           ll.plot3(lon,lat,height,varargin{:});
           ll.run;
        end
        
        function [lat,lon,height] = toWgs84Egm96(obj)
            [lat,lon,height] = ecef2geodetic(wgs84Ellipsoid,obj.x,obj.y,obj.z);
            undulationfile=util.sys.config('egm96undulation','file');
            undulation=util.io.getUndulationData(lat,lon,undulationfile);
            height = height-undulation;
        end
        
        function [lat,lon,height] = toGeodetic(obj)
            [lat,lon,height] = ecef2geodetic(wgs84Ellipsoid,obj.x,obj.y,obj.z);
        end
        
        function [e,n,u] = toEnu(obj,refpos)
            [lat,lon,height] = ecef2geodetic(wgs84Ellipsoid,refpos.x,refpos.y,refpos.z);
            [e,n,u] = ecef2enu(obj.x,obj.y,obj.z,lat,lon,height,wgs84Ellipsoid);
            if (nargout==1)
                e = util.math.Vector3(e, n, u);
            end
        end
        
        function [e,n,height] = toUtmDHHN2016(obj)
            [lat,lon,height] = ecef2geodetic(wgs84Ellipsoid,obj.x,obj.y,obj.z);
            [e,n] = utmups_fwd(lat,lon,32);
            undulationfile=util.sys.config('gcg2016undulation','file');
            undulation=util.io.getUndulationData(lat,lon,undulationfile);
            height = height-undulation;
        end
        
        function [e,n,height] = toUtm(obj)
            [lat,lon,height] = ecef2geodetic(wgs84Ellipsoid,obj.x,obj.y,obj.z);
            [e,n] = utmups_fwd(lat,lon,32);
        end
        
        function objo = toFrame(obj,frame)
            inputCoordinates = double(obj);
            outputCoordinates = util.math.transformITRF(inputCoordinates,obj.epoch,obj.referenceFrame,frame);
            objo = util.math.GeoVector(outputCoordinates);
            objo.epoch = obj.epoch;
            objo.referenceFrame = frame;
            objo.referenceEllipsoid = obj.referenceEllipsoid;
        end
        
        function objo = toEci(obj,timeSinceEpoch)
            rotationAngle =  7.2921151467e-5*timeSinceEpoch;
            objo = util.math.GeoVector(obj.x.*cos(rotationAngle)-obj.y.*sin(rotationAngle),...
                obj.x.*sin(rotationAngle)+obj.y.*cos(rotationAngle),...
                obj.z);
            objo.epoch = obj.epoch;
            objo.referenceFrame = obj.referenceFrame;
            objo.referenceEllipsoid = obj.referenceEllipsoid;
            
        end
    end
    
    methods (Static)
        function obj = fromGeodetic(lat,lon,height)
            [x,y,z] = geodetic2ecef(wgs84Ellipsoid,lat,lon,height);
            obj = util.math.GeoVector(x,y,z);
        end
        
        function obj = fromEnu(refpos,e,n,u)
            [lat,lon,height] = ecef2geodetic(wgs84Ellipsoid,refpos.x,refpos.y,refpos.z);
            [x,y,z] = enu2ecef(e,n,u,lat,lon,height,wgs84Ellipsoid);            
            obj = util.math.GeoVector(x,y,z);
        end
        
        function obj = fromWgs84Egm96(lat,lon,height)
            undulationfile=util.sys.config('egm96undulation','file');
            undulation=util.io.getUndulationData(lat,lon,undulationfile);
            [x,y,z] = geodetic2ecef(wgs84Ellipsoid,lat,lon,height+undulation);
            obj = util.math.GeoVector(x,y,z,'WGS84 (G1762)');
        end
        
        function obj = fromUtmDHHN2016(e,n,height)
            undulationfile=util.sys.config('gcg2016undulation','file');
            [lat,lon] = utmups_inv(e,n,32,1);
            undulation=util.io.getUndulationData(lat,lon,undulationfile);
            [x,y,z] = geodetic2ecef(referenceEllipsoid('Geodetic Reference System 1980'),lat,lon,height+undulation);
            obj = util.math.GeoVector(x,y,z,'ETRS89dref2016');
        end
        
        function obj = fromUtm(e,n,height)
            [lat,lon] = utmups_inv(e,n,32,1);
            [x,y,z] = geodetic2ecef(wgs84Ellipsoid, lat, lon, height);
            obj = util.math.GeoVector(x,y,z);
        end
        
    end
end

