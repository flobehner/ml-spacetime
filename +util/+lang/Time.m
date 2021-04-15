classdef Time < matlab.mixin.internal.MatrixDisplay
    %GPSTIME Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Access=private)
        datasecond int64
        datasubsecond double
    end
    
    properties (Dependent)
        utcLeapSeconds
    end
    
    methods
        function obj = Time(varargin)
            %GPSTIME Construct an instance of this class
            %   Detailed explanation goes here
            if nargin==0
                dt = datetime('now','TimeZone','UTCLeapSeconds');
                obj.datasecond = int64(floor(seconds(dt-datetime(1980,1,6,'TimeZone','UTCLeapSeconds'))));
                obj.datasubsecond = mod(dt.Second,1);
            elseif nargin==1
                if isa(varargin{1},'datetime')
                    dt = varargin{1};
                    dt.TimeZone = 'UTCLeapSeconds';
                    obj.datasecond = int64(floor(seconds(dt-datetime(1980,1,6,'TimeZone','UTCLeapSeconds'))));
                    obj.datasubsecond = mod(dt.Second,1);
                elseif isa(varargin{1},'string')
                    dt = regexp(varargin{1},"([0-9]{4})[-]?([0-9]{2})[-]?([0-9]{2})T([0-9]{2})[:]?([0-9]{2})[:]?([0-9]{2})(.[0-9])*",'tokens','once');
                    
                    if numel(varargin{1})>1
                        ds = zeros(size(varargin{1}),'int64');
                        dss = zeros(size(varargin{1}));
                        for ind=1:numel(varargin{1})
                            dpp = double(dt{ind});
                            dpp(isnan(dpp)) = 0;
                            ds(ind) = seconds(datetime(dpp(1),dpp(2),dpp(3),dpp(4),dpp(5),dpp(6),'TimeZone','UTCLeapSeconds')-datetime(1980,1,6,'TimeZone','UTCLeapSeconds'));
                            dss(ind) = dpp(7);
                            
                        end
                        obj.datasecond = ds;
                        obj.datasubsecond = dss;
                    else
                        dpp = double(dt);
                        obj.datasecond  = seconds(datetime(dpp(1),dpp(2),dpp(3),dpp(4),dpp(5),dpp(6),'TimeZone','UTCLeapSeconds')-datetime(1980,1,6,'TimeZone','UTCLeapSeconds'));
                        dpp(isnan(dpp)) = 0;
                        obj.datasubsecond = dpp(7);
                    end
                    
                else
                    obj.datasecond = int64(varargin{1});
                    obj.datasubsecond=zeros(size(varargin{1}));
                end
            else
                if any(size(varargin{1})~=size(varargin{2})) && ~isscalar(varargin{1})
                    error('Time:subsref',...
                        'Size of Second and Subseconds must be equal')
                end
                obj.datasecond = int64(varargin{1})+int64(floor(varargin{2}));
                obj.datasubsecond = mod(varargin{2},1);
            end
        end
    end
     methods (Access = protected)
        function ss = displayScalarObject(obj)
            fprintf('id: %s\n', string(obj));
        end
    end  
    methods
        
        % Display and strings
        function disp(obj)
            formatSpec = 'GPST%u:%09uns %s\n';
            stro = string(datetime(obj));
            ds = obj.datasecond;
            dss = round(obj.datasubsecond*1e9);
            ds = ds+int64(floor(dss/1e9));
            dss = mod(dss,1e9);
            fprintf("%s",strjoin(compose(formatSpec,ds(:),dss(:),stro(:))))
        end
      


        function sref = subsref(obj,s)
            % obj(i) is equivalent to obj.Data(i)
            switch s(1).type
                case '.'
                    sref = builtin('subsref',obj,s);
                case '()'
                    if length(s)<2
                        dsec = builtin('subsref',obj.datasecond,s);
                        dssec = builtin('subsref',obj.datasubsecond,s);
                        sref = util.lang.Time(dsec,dssec);
                        return
                    else
                        sref = builtin('subsref',obj,s);
                    end
                case '{}'
                    error('Time:subsref',...
                        'Not a supported subscripted reference')
            end
        end
        
        function obj = subsasgn(obj,s,val)
            if isempty(s) && isa(val,'util.lang.Time')
                obj = util.lang.Time(val.datasecond,val.datasubsecond);
            end
            switch s(1).type
                case '.'
                    obj = builtin('subsasgn',obj,s,val);
                case '()'
                    %
                    if length(s)<2
                        if isa(val,'util.lang.Time')
                            snew = substruct('.','datasecond','()',s(1).subs(:));
                            obj = subsasgn(obj,snew,val.datasecond);
                            snew = substruct('.','datasubsecond','()',s(1).subs(:));
                            obj = subsasgn(obj,snew,val.datasubsecond);
                        end
                        if isempty(val)
                            snew = substruct('.','datasecond','()',s(1).subs(:));
                           obj = subsasgn(obj,snew,[]);
                           snew = substruct('.','datasubsecond','()',s(1).subs(:));
                           obj = subsasgn(obj,snew,[]);
                        end
                    end
                case '{}'
                    error('Time:subsasgn',...
                        'Not a supported subscripted assignment')
            end
        end
        
        function ind = end(obj,k,n)
            szd = size(obj.datasecond);
            if k < n
                ind = szd(k);
            else
                ind = prod(szd(k:end));
            end
        end
        
        function ms = plus(obj,b)
            ms = util.lang.Time(obj.datasecond+int64(floor(b)),obj.datasubsecond+mod(b,1));
        end
        
        function ms = minus(obj,b)
            if isa(b,'util.lang.Time')
                dssec = obj.datasubsecond-b.datasubsecond;
                dsec = obj.datasecond-b.datasecond;
                ms = double(dsec)+dssec;
            else
                ms = util.lang.Time(obj.datasecond+int64(floor(-b)),obj.datasubsecond+mod(-b,1));
            end
        end
        
        function mv = vertcat(obj,varargin)
            ds = cellfun(@(x) x.datasecond,varargin,'UniformOutput',false);
            dss = cellfun(@(x) x.datasubsecond,varargin,'UniformOutput',false);
            dsc = vertcat(obj.datasecond,ds{:});
            dssc = vertcat(obj.datasubsecond,dss{:});
            mv = util.lang.Time(dsc,dssc);
        end
        
        function mv = horzcat(obj,varargin)
            ds = cellfun(@(x) x.datasecond,varargin,'UniformOutput',false);
            dss = cellfun(@(x) x.datasubsecond,varargin,'UniformOutput',false);
            dsc = horzcat(obj.datasecond,ds{:});
            dssc = horzcat(obj.datasubsecond,dss{:});
            mv = util.lang.Time(dsc,dssc);
        end
        
        function mv = mean(obj,varargin)
            t0 = util.lang.Time(obj.datasecond(1),obj.datasubsecond(1));
            tmp = obj-t0;
            tmp = mean(tmp,varargin{:});
            mv = t0+tmp;
        end
        
        function [mv,mi] = min(obj,varargin)
            t0 = util.lang.Time(obj.datasecond(1),obj.datasubsecond(1));
            tmp = obj-t0;
            [tmp,mi] = builtin('min',tmp,varargin{:});
            mv = t0+tmp;
        end
        
        function [mv,mi] = max(obj,varargin)
            t0 = util.lang.Time(obj.datasecond(1),obj.datasubsecond(1));
            tmp = obj-t0;
            [tmp,mi] = builtin('max',tmp,varargin{:});
            mv = t0+tmp;
        end
                
        function ms = eq(obj,b)
            ms = and(eq(obj.datasecond,b.datasecond),eq(obj.datasubsecond,b.datasubsecond));
        end
        
        function df = diff(obj)
            dssec = diff(obj.datasubsecond);
            dsec = diff(obj.datasecond);
            df = double(dsec)+dssec;
        end
        
        function sz = size(obj,varargin)
            sz = size(obj.datasecond,varargin{:});
        end
        
        function sz = datetime(obj,varargin)
            sz = (datetime(1980,1,6,'TimeZone','UTCLeapSeconds')+seconds(obj.datasecond))+seconds(obj.datasubsecond);
            sz = datetime(datetime(sz),'TimeZone','UTC');
        end
        
        function sz = posixtime(obj)
            sz = 315964800+double(obj);
        end
        
        function sz = string(obj,varargin)
            if isempty(varargin)
                varargin{1} = 'yyyyMMdd''T''HHmmss.SSSX';
            end
            sz = string(datetime(datetime(obj),'TimeZone','UTC'),varargin{:});
        end
        
        function sz = double(obj)
            sz = double(obj.datasecond)+obj.datasubsecond;
        end
        
        function sz = int64(obj)
            sz = obj.datasecond+int64(obj.datasubsecond);
        end
        
        function sz = colon(j,i,k)
            sz = j+(0:i:(k-j));
        end
        
        function sz = length(obj)
            sz = length(obj.datasecond);
        end
        
        function [sz,ind] = sort(obj,varargin)
            [sz,ind] = sort(double(obj),varargin{:});
        end
        
        function dd = get.datasubsecond(obj)
            dd = obj.datasubsecond;
        end
        
        function obj = floor(obj)
            obj = util.lang.Time(obj.datasecond);
        end
        
        function obj = transpose(obj)
            obj = util.lang.Time(obj.datasecond.',obj.datasubsecond.');
        end
        
        function obj = getStartOfGpsWeek(obj)
            gpsWeek = floor(double(obj.datasecond)/604800);
            obj = util.lang.Time(gpsWeek*604800);
        end
        
        function obj = getSecondsOfGpsWeek(obj)
            gpsWeek = floor(double(obj.datasecond)/604800);
            obj = obj-util.lang.Time(gpsWeek*604800);
        end
        
        function gpsWeek = getGpsWeek(obj)
            gpsWeek = floor(double(obj.datasecond)/604800);
        end
        
        function obj = getStartOfUtcWeek(obj)
            gpsWeek = floor(double(obj.datasecond)/604800);
            obj = util.lang.Time(gpsWeek*604800,0);
            obj = obj+obj.utcLeapSeconds;
        end
        
        function obj = getStartOfYear(obj)
            soy = dateshift(datetime(obj),'start','year');
            obj = util.lang.Time(soy);
        end
        
        function doy = getDayOfYear(obj)
            doy = (obj-getStartOfYear(obj))/86400;
        end
        
        function ls = get.utcLeapSeconds(obj)
            leapSecondTime = [46828799 78364800 109900801 173059202 252028803 315187204 346723205 393984006 425520007 457056008 504489609 551750410 599184011 820108812 914803213 1025136014 1119744015 1167264016];
            leapSecondValue = [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1];
            ls = sum(leapSecondValue.*(leapSecondTime<=obj.datasecond),2);
        end
    end
    
    methods (Static)
        function z = zeros(varargin)
            z = util.lang.Time(zeros(varargin{:}));
        end
        
        function obj = fromUnixTime(second, subsecond)
            obj = util.lang.Time(second-315964800,subsecond);
            obj = obj+obj.utcLeapSeconds;            
        end
        
        function obj = fromTerrestrialTime(mjd, timeofday)
            obj = (util.lang.Time(630763148,0.816)+86400*mjd)+timeofday;
        end
        
        function obj = fromGpsTime(week, seconds)
            obj = util.lang.Time(604800*double(week),double(seconds));
        end
        
        function obj = fromSinexEpoch(str)
          timeComponents = double(split(str,':'));
          if size(timeComponents,2) == 1
              timeComponents = timeComponents.';
          end
          ct21 = timeComponents(:,1)<50;
          timeComponents(ct21,1) = timeComponents(ct21,1)+2000;
          timeComponents(~ct21,1) = timeComponents(~ct21,1)+1900;
          timeDatetime = datetime(timeComponents(:,1),1,1)+days(timeComponents(:,2)-1)+seconds(timeComponents(:,3));
          obj = util.lang.Time(timeDatetime);
        end
        
    end
    
end

