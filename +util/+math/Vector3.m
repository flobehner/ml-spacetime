classdef Vector3
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Access=protected)
        x_
        y_
        z_
    end
    
    properties (Dependent)
        norm
        x
        y
        z
    end
    
    methods
        function obj = Vector3(varargin)
            %UNTITLED2 Construct an instance of this class
            %   Detailed explanation goes here
            if (nargin == 0)
                obj.x_ = zeros(0);
                obj.y_ = zeros(0);
                obj.z_ = zeros(0);
            elseif(nargin == 1)
                if (isa(varargin{1},'util.math.Quaternion'))
                    obj.x_ = varargin{1}.x_;
                    obj.y_ = varargin{1}.y_;
                    obj.z_ = varargin{1}.z_;
                elseif (isa(varargin{1},'util.math.Vector3'))
                    obj.x = varargin{1}.x;
                    obj.y = varargin{1}.y;
                    obj.z = varargin{1}.z;
                else
                    sz = size(varargin{1});
                    if (sz(1) == 0)
                        obj.x_ = zeros(0);
                        obj.y_ = zeros(0);
                        obj.z_ = zeros(0);
                    elseif (sz(1) == 1)
                        obj.x_=varargin{1};
                        obj.y_=varargin{1};
                        obj.z_=varargin{1};
                    elseif (sz(1) == 3)
                        data = varargin{1};
                        obj.x_ = zeros([sz(2:end),1]);
                        obj.y_ = zeros([sz(2:end),1]);
                        obj.z_ = zeros([sz(2:end),1]);
                        obj.x_(1:numel(data(1,:))) = data(1,:);
                        obj.y_(1:numel(data(2,:))) = data(2,:);
                        obj.z_(1:numel(data(3,:))) = data(3,:);
                    else
                        error('Vector3:dimensionMismatch', 'Vector3: Dimension missmatch!');
                    end
                end
            elseif(nargin == 3)
                sz = size(varargin{1});
                if (all(sz == size(varargin{2})) && all(sz == size(varargin{3})))
                    obj.x_ = varargin{1};
                    obj.y_ = varargin{2};
                    obj.z_ = varargin{3};
                else
                    error('Vector3:dimensionMismatch', 'Vector3: Dimension missmatch!');
                end
                
            end
        end
        
        function obj = gpuArray(obj)
            obj.x_ = gpuArray(obj.x_);
            obj.y_ = gpuArray(obj.y_);
            obj.z_ = gpuArray(obj.z_);
        end
        
        function obj = gather(obj)
            obj.x_ = gather(obj.x_);
            obj.y_ = gather(obj.y_);
            obj.z_ = gather(obj.z_);
        end
        
        function obj=set.x(obj,rhs)
            if (numel(rhs) == 1 && numel(obj.x_)>1)
                obj.x_(:) = rhs;
            else
                obj.x_ = rhs;
            end
            obj = obj.resize(size(obj.x_));
        end
        
        function v=get.x(obj)
            v = obj.x_;
        end
        
        
        function obj=set.y(obj,rhs)
            if (numel(rhs) == 1 && numel(obj.y_)>1)
                obj.y_(:) = rhs;
            else
                obj.y_ = rhs;
            end
            obj = obj.resize(size(obj.y_));
        end
        
        function v=get.y(obj)
            v = obj.y_;
        end
        
        function obj=set.z(obj,rhs)
            if (numel(rhs) == 1 && numel(obj.z_)>1)
                obj.z_(:) = rhs;
            else
                obj.z_ = rhs;
            end
            obj = obj.resize(size(obj.z_));
        end
        
        function v=get.z(obj)
            v = obj.z_;
        end
        
        function v=get.norm(obj)
            v = vecnorm(obj);
        end
        
        function v=vecnorm(obj)
            v = sqrt(obj.x_.^2+obj.y_.^2+obj.z_.^2);
        end
        
        function v=abs(obj)
            v = util.math.Vector3(abs(obj.x_),abs(obj.y_),abs(obj.z_));
        end
        
        function obj = setnorm(obj, nn)
            s = nn./vecnorm(obj);
            obj = util.math.Vector3(obj.x_.*s,obj.y_.*s,obj.z_.*s);
        end
        
        function obj=set.norm(obj, nn)
            obj=setnorm(obj,nn);
        end
        
        function r = plus(lhs, rhs)
            r = util.math.Vector3(lhs.x_ + rhs.x_,...
                lhs.y_ + rhs.y_,...
                lhs.z_ + rhs.z_);
        end
        
        function r = minus(lhs, rhs)
            r = util.math.Vector3(lhs.x_ - rhs.x_,...
                lhs.y_ - rhs.y_,...
                lhs.z_ - rhs.z_);
        end
        
        function r = uminus(rhs)
            r = util.math.Vector3(- rhs.x_,...
                -rhs.y_,...
                -rhs.z_);
        end
        
        function obj = transpose(obj)
            obj = feval(class(obj),obj.x_.',obj.y_.',obj.z_.');
        end
        
        function obj = reshape(obj,varargin)
            obj = util.math.Vector3(reshape(obj.x_,varargin{:}),reshape(obj.y_,varargin{:}),reshape(obj.z_,varargin{:}));
        end
        
        function obj = permute(obj,order)
            obj = util.math.Vector3(permute(obj.x_,order),permute(obj.y_,order),permute(obj.z_,order));
        end
        
        function r = times(lhs, rhs)
            % if (isscalar(lhs) || isscalar(rhs) || (all(size(lhs) == size(rhs))))
            if (isa(lhs,'util.math.Vector3') & isa(rhs,'util.math.Vector3'))
                r = util.math.Vector3(lhs.x_.*rhs.x_, lhs.y_.*rhs.y_, lhs.z_.*rhs.z_);
            elseif (~isa(lhs,'util.math.Vector3'))
                r = util.math.Vector3(lhs.*rhs.x_, lhs.*rhs.y_, lhs.*rhs.z_);
            elseif (~isa(rhs,'util.math.Vector3'))
                r = util.math.Vector3(lhs.x_.*rhs, lhs.y_.*rhs, lhs.z_.*rhs);
            end
            %r = util.math.Vector3(double(lhs) .* double(rhs));
            %else
            %    r = util.math.Vector3(double(lhs).*double(rhs));
            %end
        end
        
        function r = rdivide(lhs, rhs)
            r = util.math.Vector3(lhs.x_./rhs,lhs.y_./rhs,lhs.z_./rhs);
        end
        
        function r = cross(lhs, rhs)
            r = util.math.Vector3(lhs.y_.*rhs.z_-lhs.z_.*rhs.y_,...
                lhs.z_.*rhs.x_-lhs.x_.*rhs.z_,...
                lhs.x_.*rhs.y_-lhs.y_.*rhs.x_);
        end
        
        function r = dot(lhs,rhs)
            r = lhs.x_.*rhs.x_ + lhs.y_.*rhs.y_ + lhs.z_.*rhs.z_;
        end
        
        function r = mtimes(lhs, rhs)
            if (isscalar(lhs) || isscalar(rhs) || (all(size(lhs) == size(rhs)) && (~isa(lhs,'util.math.Vector3') || ~isa(rhs,'util.math.Vector3'))))
                r = lhs.*rhs;
            elseif (isa(rhs,'util.math.Quaternion'))
                r = rhs*lhs;
            elseif (isa(lhs,'util.math.Vector3') && isa(rhs,'util.math.Vector3'))
                r = dot(lhs,rhs);
            elseif (isa(rhs,'util.math.Vector3'))
                r = util.math.Vector3(double(lhs)*rhs.data);
            else
                error('Vector3:dimensionMismatch', 'Vector3: Dimension missmatch!');
            end
            
        end
        
        function r = isscalar(obj)
            r = false;
        end
        
        function ind = end(obj,k,n)
            szd = size(obj);
            if k < n
                ind = szd(k);
            else
                ind = prod(szd(k:end));
            end
        end
        
        function d = double(obj)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            d = cat(1,shiftdim(obj.x_,-1),shiftdim(obj.y_,-1),shiftdim(obj.z_,-1));
        end
        
        function n = numArgumentsFromSubscript(obj,s,indexingContext)
            switch s(1).type
                case '()'
                    n = 1;
                otherwise
                    n = builtin('numArgumentsFromSubscript', obj, s, indexingContext);
            end
        end
        
        function obj  = subsref(obj,s)
            % obj(i) is equivalent to obj.Data(i)
            switch s(1).type
                case '.'
                    obj = builtin('subsref',obj,s);
                case '()'
                    obj = util.math.Vector3(builtin('subsref',obj.x_,s(1)),builtin('subsref',obj.y_,s(1)),builtin('subsref',obj.z_,s(1)));
                    if (length(s) > 1)
                        obj = obj.subsref(s(2:end));
                    end
                case '{}'
                    error('Vector3:subsref',...
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
                            error('Vector3:subsasgn',...
                                'Not a supported subscripted reference');
                        end
                    else
                        a = util.math.Vector3(a);
                        obj = util.math.Vector3(obj);
                        obj.x = builtin('subsasgn',obj.x_,s,a.x_);
                        obj.y = builtin('subsasgn',obj.y_,s,a.y_);
                        obj.z = builtin('subsasgn',obj.z_,s,a.z_);
                    end
                case '{}'
                    error('Vector3:subsasgn',...
                        'Not a supported subscripted reference')
            end
        end
        
        function disp(obj,name)
            if (nargin==1)
                name = '';
            end
            disp(double(obj));
        end
        
        function s=size(obj,varargin)
            if builtin('length',obj)==1
                s=size(obj.x_,varargin{:});
            else
                s = arrayfun(@(obj) size(obj.x_,varargin{:}),obj);
            end
        end
        
        function l=length(obj)
            s=size(obj.x_);
            l = max(s);
        end
        
        function l=numel(obj)
            l = numel(obj.x);
        end
        
        function ph = plot3(obj,varargin)
            ph=plot3(obj.x_,obj.y_,obj.z_,varargin{:});
        end
        
        function mv = vertcat(obj,varargin)
            xc = cellfun(@(x) x.x_,varargin,'UniformOutput',false);
            yc = cellfun(@(x) x.y_,varargin,'UniformOutput',false);
            zc = cellfun(@(x) x.z_,varargin,'UniformOutput',false);
            x = vertcat(obj.x_,xc{:});
            y = vertcat(obj.y_,yc{:});
            z = vertcat(obj.z_,zc{:});
            mv = util.math.Vector3(x,y,z);
        end
        
        function mv = horzcat(obj,varargin)
            xc = cellfun(@(x) x.x_,varargin,'UniformOutput',false);
            yc = cellfun(@(x) x.y_,varargin,'UniformOutput',false);
            zc = cellfun(@(x) x.z_,varargin,'UniformOutput',false);
            x = horzcat(obj.x_,xc{:});
            y = horzcat(obj.y_,yc{:});
            z = horzcat(obj.z_,zc{:});
            mv = util.math.Vector3(x,y,z);
        end
        
        
        function mv = sum(obj,varargin)
            x= sum(obj.x_,varargin{:});
            y= sum(obj.y_,varargin{:});
            z= sum(obj.z_,varargin{:});
            mv = util.math.Vector3(x,y,z);
        end
        
        function mv = diff(obj,varargin)
            x= diff(obj.x_,varargin{:});
            y= diff(obj.y_,varargin{:});
            z= diff(obj.z_,varargin{:});
            mv = util.math.Vector3(x,y,z);
        end
        
        function mv = gradient(obj,varargin)
            x= gradient(obj.x_,varargin{:});
            y= gradient(obj.y_,varargin{:});
            z= gradient(obj.z_,varargin{:});
            mv = util.math.Vector3(x,y,z);
        end
        
        function mv = cumtrapz(obj,varargin)
            if length(varargin)>=1
                x= cumtrapz(varargin{1},obj.x_,varargin{2:end});
                y= cumtrapz(varargin{1},obj.y_,varargin{2:end});
                z= cumtrapz(varargin{1},obj.z_,varargin{2:end});
            else
                x= cumtrapz(obj.x_);
                y= cumtrapz(obj.y_);
                z= cumtrapz(obj.z_);
            end
            mv = util.math.Vector3(x,y,z);
        end
        
    end
    
    methods(Access=protected)
        
        function obj = resize(obj,sz)
            for prop = ["x_","y_","z_"]
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
    
    methods (Static)
        function obj = fromSpherical(theta,phi,r)
            rst = r.*sin(theta);
           obj = util.math.Vector3(cos(phi).*rst,...
               sin(phi).*rst,...
               r.*cos(theta));
            
        end
    end
end

