classdef Quaternion < util.math.Vector3
    %UNTITLED4 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(Access=protected)
        w_
    end
    
    properties(Dependent)
        w
    end
    
    
    methods
        function obj = Quaternion(varargin)
            obj@util.math.Vector3();
            if (nargin == 0)
                obj.w_ = zeros(0);
                obj.x_ = zeros(0);
                obj.y_ = zeros(0);
                obj.z_ = zeros(0);
            elseif(nargin == 1)
                if (isa(varargin{1},'util.math.Quaternion'))
                    obj = varargin{1};
                elseif (isa(varargin{1},'util.math.Vector3'))
                    obj.x_ = varargin{1}.x_;
                    obj.y_ = varargin{1}.y_;
                    obj.z_ = varargin{1}.z_;
                    obj.w_ = zeros(size(obj.x_));
                else
                    sz = size(varargin{1});
                    if (sz(1) == 0)
                        obj.w_ = zeros(0);
                        obj.x_ = zeros(0);
                        obj.y_ = zeros(0);
                        obj.z_ = zeros(0);
                    elseif (sz(1) == 1)
                        obj.w_=varargin{1};
                        obj.x_=varargin{1};
                        obj.y_=varargin{1};
                        obj.z_=varargin{1};
                    elseif (sz(1) == 3)
                        data = varargin{1};
                        obj.w_ = zeros([sz(2:end),1]);
                        obj.x_ = zeros([sz(2:end),1]);
                        obj.y_ = zeros([sz(2:end),1]);
                        obj.z_ = zeros([sz(2:end),1]);
                        obj.w_(1:numel(data(1,:))) = 0;
                        obj.x_(1:numel(data(2,:))) = data(1,:);
                        obj.y_(1:numel(data(3,:))) = data(2,:);
                        obj.z_(1:numel(data(4,:))) = data(3,:);
                    elseif (sz(1) == 4)
                        data = varargin{1};
                        obj.w_ = zeros([sz(2:end),1]);
                        obj.x_ = zeros([sz(2:end),1]);
                        obj.y_ = zeros([sz(2:end),1]);
                        obj.z_ = zeros([sz(2:end),1]);
                        obj.w_(1:numel(data(1,:))) = data(1,:);
                        obj.x_(1:numel(data(2,:))) = data(2,:);
                        obj.y_(1:numel(data(3,:))) = data(3,:);
                        obj.z_(1:numel(data(4,:))) = data(4,:);
                        
                    else
                        error('Vector3:dimensionMismatch', 'Vector3: Dimension missmatch!');
                    end
                end
            elseif(nargin == 2)
                if (isa(varargin{1},'util.math.Vector3'))
                    obj.x_ = varargin{1}.x_;
                    obj.y_ = varargin{1}.y_;
                    obj.z_ = varargin{1}.z_;
                    obj.w_ = varargin{2};
                end
            elseif(nargin == 3)
                sz = size(varargin{1});
                if (all(sz == size(varargin{2})) && all(sz == size(varargin{3})))
                    obj.x_ = varargin{1};
                    obj.y_ = varargin{2};
                    obj.z_ = varargin{3};
                    obj.w_ = zeros(size(obj.x_));
                else
                    error('Vector3:dimensionMismatch', 'Vector3: Dimension missmatch!');
                end
            elseif(nargin == 4)
                sz = size(varargin{1});
                if (all(sz == size(varargin{2})) && all(sz == size(varargin{3})) && all(sz == size(varargin{4})))
                    obj.w_ = varargin{1};
                    obj.x_ = varargin{2};
                    obj.y_ = varargin{3};
                    obj.z_ = varargin{4};
                else
                    error('Vector3:dimensionMismatch', 'Vector3: Dimension missmatch!');
                end
            end
        end
        
        function obj = gpuArray(obj)
            obj.x_ = gpuArray((obj.x_));
            obj.y_ = gpuArray((obj.y_));
            obj.z_ = gpuArray((obj.z_));
            obj.w_ = gpuArray((obj.w_));
        end
        
        function obj = gather(obj)
            obj.x_ = gather(obj.x_);
            obj.y_ = gather(obj.y_);
            obj.z_ = gather(obj.z_);
            obj.w_ = gather(obj.w_);
        end
        
        function obj=set.w(obj,rhs)
            if (numel(rhs) == 1 && numel(obj.w_)>1)
                obj.w_(:) = rhs;
            else
                obj.w_ = rhs;
            end
            obj = obj.resize(size(obj.w_));
        end
        
        function v=get.w(obj)
            v = obj.w_;
        end
        
        function v=vecnorm(obj)
            v = sqrt(obj.x_.^2+obj.y_.^2+obj.z_.^2+obj.w_.^2);
        end
        
        function obj = setnorm(obj, nn)
            s = nn./vecnorm(obj);
            obj = util.math.Quaternion(obj.w_.*s,obj.x_.*s,obj.y_.*s,obj.z_.*s);
        end
        
        function r = plus(lhs, rhs)
            r = util.math.Vector3(...
                lhs.w_ + rhs.w_,...
                lhs.x_ + rhs.x_,...
                lhs.y_ + rhs.y_,...
                lhs.z_ + rhs.z_);
        end
        
        function r = minus(lhs, rhs)
            r = util.math.Vector3(...
                lhs.w_ - rhs.w_,...
                lhs.x_ - rhs.x_,...
                lhs.y_ - rhs.y_,...
                lhs.z_ - rhs.z_);
        end
        
        
        function r = dot(lhs,rhs)
            r = lhs.x_.*rhs.x_ + lhs.y_.*rhs.y_ + lhs.z_.*rhs.z_ + lhs.w_.*rhs.w_;
        end
        
        function obj = conj(obj)
            obj.x_ = -obj.x_;
            obj.y_ = -obj.y_;
            obj.z_ = -obj.z_;
        end
        
        function r = times(lhs, rhs)
            if (isa(lhs,'util.math.Quaternion') & isa(rhs,'util.math.Quaternion'))
                r = util.math.Vector3(lhs.w_.*rhs.w_, lhs.x_.*rhs.x_, lhs.y_.*rhs.y_, lhs.z_.*rhs.z_);
            elseif (~isa(lhs,'util.math.Quaternion'))
                r = util.math.Vector3(lhs.*rhs.w_, lhs.*rhs.x_, lhs.*rhs.y_, lhs.*rhs.z_);
            elseif (~isa(rhs,'util.math.Quaternion'))
                r = util.math.Vector3(lhs.w_.*rhs, lhs.x_.*rhs, lhs.y_.*rhs, lhs.z_.*rhs);
            end
            %             if (isscalar(lhs) || isscalar(rhs) || (all(size(lhs) == size(rhs)) && isa(lhs,'util.math.Quaternion') && isa(rhs,'util.math.Quaternion')))
            %                 r = util.math.Quaternion(double(lhs) .* double(rhs));
            %             else
            %                 if (all(size(rhs) == size(lhs)))
            %                     if (isa(rhs,'util.math.Quaternion'))
            %                         r = util.math.Quaternion(shiftdim(double(lhs),-1).* double(rhs));
            %                     else
            %                         r = util.math.Quaternion(double(lhs) .* shiftdim(double(rhs),-1));
            %                     end
            %
            %                 else
            %                     r = util.math.Quaternion(double(lhs).*double(rhs));
            %                 end
            %             end
        end
        
        function obj = inv(obj)
            obj = conj(obj)./obj.norm.^2;
        end
        
        function r = rdivide(lhs, rhs)
            r = util.math.Quaternion(lhs.w_./rhs,lhs.x_./rhs,lhs.y_./rhs,lhs.z_./rhs);
        end
        
        function [roll, pitch, yaw] = toEuler(obj)
            roll = atan2d(2*(obj.w.*obj.x+obj.y.*obj.z),1-2*(obj.x.^2+obj.y.^2));
            pitch = asind(2*(obj.w.*obj.y-obj.z.*obj.x));
            yaw = atan2d(2*(obj.w.*obj.z+obj.x.*obj.y),1-2*(obj.y.^2+obj.z.^2));
        end
        
        function obj = transpose(obj)
            obj = util.math.Quaternion(obj.w_.',obj.x_.',obj.y_.',obj.z_.');
        end
        
        function obj = reshape(obj,varargin)
            obj = util.math.Quaternion(reshape(obj.w_,varargin{:}),reshape(obj.x_,varargin{:}),reshape(obj.y_,varargin{:}),reshape(obj.z_,varargin{:}));
        end
        
        function obj = permute(obj,order)
            obj = util.math.Quaternion(permute(obj.w_,order),permute(obj.x_,order),permute(obj.y_,order),permute(obj.z_,order));
        end
        
        function r = mtimes(lhs,rhs)
            if (isa(lhs,'util.math.Vector3') && ~isa(lhs,'util.math.Quaternion')) % => Rotation
                x0 = rhs.w_; x1 = rhs.x_; x2 = rhs.y_; x3 = rhs.z_;
                x00 = x0.^2; x11 = x1.^2; x22 = x2.^2; x33 = x3.^2;
                x01 = x0.*x1; x02 = x0.*x2; x03 = x0.*x3;
                x12 = x1.*x2; x13 = x1.*x3;
                x23 = x2.*x3;
                
                y1 = lhs.x_; y2 = lhs.y_; y3 = lhs.z_;
                r = util.math.Vector3(...
                    2*(x02.*y3 + x13.*y3 - x03.*y2 + x12.*y2) + y1.*(x00 + x11 - x22 - x33),...
                    2*(x03.*y1 + x12.*y1 - x01.*y3 + x23.*y3) + y2.*(x00 - x11 + x22 - x33),...
                    2*(x01.*y2 - x02.*y1 + x13.*y1 + x23.*y2) + y3.*(x00 - x11 - x22 + x33)...
                    );
                return
            end
            
            if (isa(rhs,'util.math.Vector3') && ~isa(rhs,'util.math.Quaternion')) % => Rotation
                x0 = lhs.w_; x1 = lhs.x_; x2 = lhs.y_; x3 = lhs.z_;
                x00 = x0.^2; x11 = x1.^2; x22 = x2.^2; x33 = x3.^2;
                x01 = x0.*x1; x02 = x0.*x2; x03 = x0.*x3;
                x12 = x1.*x2; x13 = x1.*x3;
                x23 = x2.*x3;
                
                y1 = rhs.x_; y2 = rhs.y_; y3 = rhs.z_;
                r = util.math.Vector3(...
                    2*(x02.*y3 + x13.*y3 - x03.*y2 + x12.*y2) + y1.*(x00 + x11 - x22 - x33),...
                    2*(x03.*y1 + x12.*y1 - x01.*y3 + x23.*y3) + y2.*(x00 - x11 + x22 - x33),...
                    2*(x01.*y2 - x02.*y1 + x13.*y1 + x23.*y2) + y3.*(x00 - x11 - x22 + x33)...
                    );
                
                return
            end
            
            
            if (~isa(lhs,'util.math.Quaternion'))
                lhs = util.math.Quaternion(lhs);
            end
            if (~isa(rhs,'util.math.Quaternion'))
                rhs = util.math.Quaternion(rhs);
            end
            if (isa(lhs,'util.math.Quaternion') && isa(rhs,'util.math.Quaternion'))
                x0 = lhs.w_; x1 = lhs.x_; x2 = lhs.y_; x3 = lhs.z_;
                y0 = rhs.w_; y1 = rhs.x_; y2 = rhs.y_; y3 = rhs.z_;
                r = util.math.Quaternion(...
                    x0.*y0-x1.*y1-x2.*y2-x3.*y3,...
                    x0.*y1+x1.*y0+x2.*y3-x3.*y2,...
                    x0.*y2-x1.*y3+x2.*y0+x3.*y1,...
                    x0.*y3+x1.*y2-x2.*y1+x3.*y0);
            else
                error('Quaternion:typeMissmatch', 'Quaternion: Type missmatch!');
            end
        end
        
        function n = numArgumentsFromSubscript(obj,s,indexingContext)
            n = 1;
            switch s(1).type
                case '()'
                    n = 1;
                otherwise
                    n = builtin('numArgumentsFromSubscript', obj, s, indexingContext);
            end
        end
        
        function sref = subsref(obj,s)
            % obj(i) is equivalent to obj.Data(i)
            switch s(1).type
                case '.'
                    sref = builtin('subsref',obj,s);
                case '()'
                    sref = util.math.Quaternion(builtin('subsref',obj.w_,s(1)),...
                        builtin('subsref',obj.x_,s(1)),builtin('subsref',obj.y_,s(1)),builtin('subsref',obj.z_,s(1)));
                    if (length(s) > 1)
                        sref = sref.subsref(s(2:end));
                    end
                case '{}'
                    error('Quaternion:subsref',...
                        'Not a supported subscripted reference')
            end
        end
        
        function obj = subsasgn(obj,s, a)
            % obj(i) is equivalent to obj.Data(i)
            switch s(1).type
                case '.'
                    obj = builtin('subsasgn',obj,s,a);
                case '()'
                    if (length(s) > 1)
                        if (s(2).type == '.')
                            s2(1) = s(2);
                            s2(2) = s(1);
                            obj = builtin('subsasgn',obj,s2,a);
                        else
                            error('Quaternion:subsasgn',...
                                'Not a supported subscripted reference');
                        end
                    else
                        a = util.math.Quaternion(a);
                        obj = util.math.Quaternion(obj);
                        obj.w_ = builtin('subsasgn',obj.w_,s,a.w_);
                        obj.x_ = builtin('subsasgn',obj.x_,s,a.x_);
                        obj.z_ = builtin('subsasgn',obj.z_,s,a.z_);
                        obj.y_ = builtin('subsasgn',obj.y_,s,a.y_);
                    end
                case '{}'
                    error('Quaternion:subsasgn',...
                        'Not a supported subscripted reference')
            end
        end
        
        function d = double(obj)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            d = cat(1,shiftdim(obj.w_,-1),shiftdim(obj.x_,-1),shiftdim(obj.y_,-1),shiftdim(obj.z_,-1));
        end
        
        function kmlfile = showInGoogleEarth(obj,position,filename)
            if nargin<3
                filename = "temp";
            end
            kmlfile = kml([char(filename) '.kml']);
            s.type = '()';
            if size(obj)==size(position)
                for ind = 1:numel(obj)
                    s.subs = {ind};
                    kmlAddCS(kmlfile,subsref(obj,s),position(ind),'scale',0.1,'name',['Quaternion(' num2str(ind) ')']);
                end
            elseif length(obj)==1
                for ind = 1:numel(obj)
                    kmlAddCS(kmlfile,obj,position(ind));
                end
                
            elseif  length(position)==1
                for ind = 1:numel(obj)
                    kmlAddCS(kmlfile,obj(ind),position);
                end
            end
            kmlfile.run();
        end
        
    end
    
    methods(Access=protected)
        
        function obj = resize(obj,sz)
            for prop = ["x_","y_","z_","w_"]
                temp1 = obj.(prop);
                if (length(sz) ~= length(size(temp1)))
                    temp2 = zeros(sz);
                    temp2(1:numel(temp1)) = temp1(:);
                    obj.(prop) = temp2;
                elseif (sz ~= size(temp1))
                    temp2 = zeros(sz);
                    temp2(1:numel(temp1)) = temp1(:);
                    obj.(prop) = temp2;
                end
            end
        end
    end
    
    methods(Static)
        function q = fromEuler(yaw, pitch, roll)
            cy = cosd(yaw * 0.5);
            sy = sind(yaw * 0.5);
            cp = cosd(pitch * 0.5);
            sp = sind(pitch * 0.5);
            cr = cosd(roll * 0.5);
            sr = sind(roll * 0.5);
            
            q=util.math.Quaternion(cy .* cp .* cr + sy .* sp .* sr,...
                cy .* cp .* sr - sy .* sp .* cr,...
                sy .* cp .* sr + cy .* sp .* cr,...
                sy .* cp .* cr - cy .* sp .* sr);
        end
        
        function q = fromNed2Ecef(lat,lon)
            if class(lat) == "util.math.GeoVector"
                [lat,lon,~] = toGeodetic(lat);
            end
            cy = cosd(lon * 0.5);
            sy = sind(lon * 0.5);
            cp = cosd((90+lat) * 0.5);
            sp = -sind((90+lat) * 0.5);
            q=util.math.Quaternion(cy .* cp, - sy .* sp,...
                cy .* sp, sy .* cp);
        end
        
        function q = fromRotation(axis,angle)
            cy = cosd(angle * 0.5);
            sy = sind(angle * 0.5);
            q=util.math.Quaternion(axis.*sy,cy);
        end
        
        function q = fromDCM(varargin)
            if (nargin == 1)
                dcm = varargin{1};
            end
            
            if (nargin == 3)
                ex = varargin{1};
                ey = varargin{2};
                ez = varargin{3};
                
                dcm = [shiftdim(ex.x,-2), shiftdim(ey.x,-2), shiftdim(ez.x,-2);
                    shiftdim(ex.y,-2), shiftdim(ey.y,-2), shiftdim(ez.y,-2);
                    shiftdim(ex.z,-2), shiftdim(ey.z,-2), shiftdim(ez.z,-2);];
            end
            
            w = 0.5*sqrt(dcm(1,1,:)+dcm(2,2,:)+dcm(3,3,:)+1);
            x = (dcm(3,2,:)-dcm(2,3,:))./w/4;
            y = (dcm(1,3,:)-dcm(3,1,:))./w/4;
            z = (dcm(2,1,:)-dcm(1,2,:))./w/4;
            
            if(nargin==3)
                w = shiftdim(w,1);
                x = shiftdim(x,1);
                y = shiftdim(y,1);
                z = shiftdim(z,1);
            end
            q=util.math.Quaternion(w,x,y,z);
            q.norm = 1;
        end
        
        
    end
end

function target = kmlAddCS(this,attitude,position,varargin)

target = struct('type','','id','','location_id','','orientation_id','','scale_id','','model_id','');

%[longDEG,latDEG] = this.checkUnit(long,lat);

p = inputParser;

p.addParamValue('scale',1,@(a)isnumeric(a) && numel(a)==1);
p.addParamValue('color','FFFFFFFF',@(a)ischar(a) && numel(a)==8);
p.addParamValue('name','kml_quiver3D',@ischar);
p.addParamValue('id',kml.getTempID('kml_cs'),@ischar);
p.addParamValue('description','',@ischar);
p.addParamValue('visibility',true,@islogical);
p.addParamValue('model','',@ischar);
p.addParamValue('altitudeMode','relativeToGround',@(a)ismember(a,{'clampToGround','relativeToGround','absolute'}));

p.addParamValue('timeStamp','',@ischar);
p.addParamValue('timeSpanBegin','',@ischar);
p.addParamValue('timeSpanEnd','',@ischar);

p.parse(varargin{:});

arg = p.Results;


path=fileparts(mfilename('fullpath'));
arrowfileRed =['arrow3dred.dae'];
arrowfileGreen =['arrow3dgreen.dae'];
arrowfileBlue =['arrow3dblue.dae'];


modelpath = which(arrowfileRed);
if ~isempty(modelpath)
    if isempty(dir(arrowfileRed))
        copyfile(modelpath,arrowfileRed);
    end
else
    error('File %s not found!',arrowfileRed);
end
modelpath = which(arrowfileGreen);
if ~isempty(modelpath)
    if isempty(dir(arrowfileGreen))
        copyfile(modelpath,arrowfileGreen);
    end
else
    error('File %s not found!',arrowfileGreen);
end
modelpath = which(arrowfileBlue);
if ~isempty(modelpath)
    if isempty(dir(arrowfileBlue))
        copyfile(modelpath,arrowfileBlue);
    end
else
    error('File %s not found!',arrowfileBlue);
end

[lat,lon,h] = toWgs84Egm96(util.math.GeoVector(position));
nx = attitude*util.math.Vector3(1,0,0);
ny = attitude*util.math.Vector3(0,1,0);
nz = attitude*util.math.Vector3(0,0,1);
[nxE,nxN,nxU] = ecef2enuv(nx.x,nx.y,nx.z,lat,lon);
[nyE,nyN,nyU] = ecef2enuv(ny.x,ny.y,ny.z,lat,lon);
[nzE,nzN,nzU] = ecef2enuv(nz.x,nz.y,nz.z,lat,lon);

f = this.createFolder(arg.name);

if strcmpi(this.unit,'rad')
    error("wrong Conf rad");
end

headingx =  atan2d(nxE,nxN);
tiltx    =-  (atan2d(nxU,sqrt(nxE^2+nxN^2)));
rollx    = 0;
headingy =   atan2d(nyE,nyN);
tilty    =-  (atan2d(nyU,sqrt(nyE^2+nyN^2)));
rolly    = 0;
headingz =  atan2d(nzE,nzN);
tiltz    =-  (atan2d(nzU,sqrt(nzE^2+nzN^2)));
rollz    = 0;

scale = arg.scale;

target(1) = f.model(lon,lat,h-2,headingx,tiltx,rollx, 'scale',scale,'model',arrowfileRed, ...
    'altitudeMode','absolute', ...
    'name','x' ...
    );
target(2) = f.model(lon,lat,h-2,headingy,tilty,rolly, 'scale',scale,'model',arrowfileGreen, ...
    'altitudeMode','absolute', ...
    'name','y' ...
    );
target(3) = f.model(lon,lat,h-2,headingz,tiltz,rollz, 'scale',scale,'model',arrowfileBlue, ...
    'altitudeMode','absolute', ...
    'name','z' ...
    );

end



