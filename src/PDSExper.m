classdef PDSExper3 < handle
    properties % Region 0
        % Region 0 properties, default values are for FC-72
        % Thermal conductivity
        Region0_k0(1,1) double {mustBeFinite,mustBePositive} = 0.057 ;
        
        % Density
        Region0_rho0(1,1) double {mustBeFinite,mustBePositive} = 1680 ;
        
        % specific heat capacity
        Region0_c0(1,1) double {mustBeFinite,mustBePositive} = 1100 ;
        
        % Reflection coefficent
        Region0_R1(1,1) double {mustBeFinite,mustBeNonnegative,mustBeLessThanOrEqual(Region0_R1,1)} = 0 ; 
    end
    properties % Region 1
        % Region 1, Thin Film
        % Default values are for gold
        % Thermal conductivity of thin film
        Region1_k1(1,1) double {mustBeFinite,mustBePositive} = 314 ;
        
        % Density
        Region1_rho1(1,1) double {mustBeFinite,mustBePositive} = 19300 ;
        
        % Specific heat capacity
        Region1_c1(1,1) double {mustBeFinite,mustBePositive} = 129 ;
        
        % Optical absorbtion coefficent at excitation frequency
        Region1_a1(1,1) double {mustBeFinite,mustBeNonnegative} = 5*10^7 ;
        
        % Thin film thickness
        Region1_L1(1,1) double {mustBeFinite,mustBePositive} = 200*10^-9 ; 
        
        % Reflection coeff between 1 & 2
        Region1_R2(1,1) double {mustBeFinite,mustBeNonnegative,mustBeLessThanOrEqual(Region1_R2,1)} = 0 ;
        
        
    end
    properties % Region 2
        % Region 2, Substrate
        % Default values are for glass
        Region2_k2(1,1) double {mustBeFinite,mustBePositive} = 1.0 ;
        Region2_rho2(1,1) double {mustBeFinite,mustBePositive} = 2400 ;
        Region2_c2(1,1) double {mustBeFinite,mustBePositive} = 670 ;
        Region2_a2(1,1) double {mustBeFinite,mustBeNonnegative} = 0 ;
        Region2_L2(1,1) double {mustBeFinite,mustBePositive} = 1*10^-3 ;
        Region2_Rth(1,1) double {mustBeFinite,mustBeNonnegative} = 0 ; % Interfacial Thermal resistance
        
        
    end 
    properties % Region 3
        % Region 3 properties, default values are for FC-72
        % Thermal conductivity
        Region3_k3(1,1) double {mustBeFinite,mustBePositive} = 0.057 ;
        
        % Density
        Region3_rho3(1,1) double {mustBeFinite,mustBePositive} = 1680 ;
        
        % specific heat capacity
        Region3_c3(1,1) double {mustBeFinite,mustBePositive} = 1100 ;
        
        % refractive index of FC-72
        n = 1.251 ;
        
        % Thermo-optical coefficent of FC-72
        dNdT = -3.4e-4 ;
    end
    properties % Other experimental variables
        
        % Pump beam 1/e^2 radius.
        PumpRadius(1,1) double {mustBeFinite,mustBePositive} = 100e-6 ; 
        
        % Pump beam power, in W/m^2
        Power(1,1) double {mustBeFinite,mustBePositive} = 270e-06 ;
        
        % Typically the beam is not centered at zero. This shifts the
        % position so the simulation matches with the experimental data
        xshift(1,1) double {mustBeFinite} = 8.49e-3 ;
        
        % Vertically shift the model. Should always be zero.
        yshift(1,1) double {mustBeFinite} = 0 ;
        
        % FWHM of the probe beam in meters
        ProbeRadius(1,1) double {mustBeFinite,mustBePositive} = 110e-6 ;
        
        % Height above thin film(negative) or below substrate (positive)
        z(1,1) double {mustBeFinite} = -50e-6 ;
        
        % Default modulation Frequency
        f(1,1) double {mustBeFinite,mustBePositive} = 10 ;
        
        % Add an arbitrary phase shift to the model
        % this should be zero most of the time
        PhaseShift(1,1) double {mustBeFinite} = 0 ; 
        
        % Positions to calculate the model
        xdata
        
        % Power
        P = 1
    end
    properties (Access = public) % Internal variables
        %Coefficents
        EE
        B0
        Lth0
        
        M
        Lth3
        B3
        % Parameters to define how numerical integration is performed
        del_min = 0 ;
        del_max = 200000 ;
        n_regions= 5000 ;
        n_points
        del
        
        % Defines the ray mesh. Use a GenerateMeshes method to generate.
        x_raySample
        z_raySample
        x
        y
        xwidth
        zwidth
        xstep
        Weights
        LB
        UB

        % Internal variables for the precalculation of deflection values
        r_precalc
        T0dzOut_precalc
        T3dzOut_precalc
        
        % Normal deflection, without weighting
        PsiNormal
        PsiFinal
    end
    methods (Access = private) % Define internal methods
        function region0_CalculateCoeff(obj)
            % This function calculates the coefficents nessicary to compute the temeperature
            % profile of a bi-layer sample under periodic gaussian illumination.
            % 
            % More information can be found in this paper,
            % 
            % M. Reichling and H. Grönbeck, “Harmonic heat flow in isotropic layered 
            % systems and its use for thin film thermal conductivity measurements,” 
            % Journal of Applied Physics, vol. 75, no. 4, pp. 1914–1922, Feb. 1994.

            % Region 0, Fluid

            k0   = obj.Region0_k0 ;
            rho0 = obj.Region0_rho0 ;
            c0   = obj.Region0_c0 ;
            R1   = obj.Region0_R1 ;

            % Region 1, Thin Film

            k1   = obj.Region1_k1 ;
            rho1 = obj.Region1_rho1 ;
            c1   = obj.Region1_c1 ;
            a1   = obj.Region1_a1 ;
            L1   = obj.Region1_L1 ;
            R2   = obj.Region1_R2 ;

            % Region 2, Substrate

            k2   = obj.Region2_k2 ;
            rho2 = obj.Region2_rho2 ;
            c2   = obj.Region2_c2 ;
            a2   = obj.Region2_a2 ;
            L2   = obj.Region2_L2 ;
            Rth  = obj.Region2_Rth ;

            % Region 3, Fluid

            k3   = obj.Region3_k3 ;
            rho3 = obj.Region3_rho3 ;
            c3   = obj.Region3_c3 ;

            f = obj.f ; %#ok<*PROP>
            del = obj.del ;

            % Calculate
            Lth0 = sqrt( k0./(pi .* f .* rho0 .* c0)) ;
            Lth1 = sqrt( k1./(pi .* f .* rho1 .* c1)) ;
            Lth2 = sqrt( k2./(pi .* f .* rho2 .* c2)) ;
            Lth3 = sqrt( k3./(pi .* f .* rho3 .* c3)) ;

            B0 = sqrt( del.^2 + 2i./Lth0.^2 ) ;
            B1 = sqrt( del.^2 + 2i./Lth1.^2 ) ;
            B2 = sqrt( del.^2 + 2i./Lth2.^2 ) ;
            B3 = sqrt( del.^2 + 2i./Lth3.^2 ) ;

            g = k0.*B0./(k1.*B1) ;
            b = k2.*B2./(k1.*B1) ;
            v = a1./B1 ;
            n = k1.*B1.*Rth ;
            s = a2./B2 ;
            d = k3.*B3./(k2.*B2) ;
            mu = ((1-d)./(1+d)) .* exp(-2.*B2.*L2) ;
            vv = ((s-d)./(1+d)) .* exp(-(a2+B2).*L2) ;
            gamma = (1-mu)./(1+mu) ;



            H = (g+1).*(1+b.*gamma.*(1+n))-(1-b.*gamma.*(1-n)).*(1-g).*exp(-2.*B1.*L1) ;

            Gamma = (1-R1).*(a1.*obj.P).*exp(-(del.*obj.PumpRadius).^2./8)./(2.*pi.*k1)./(B1.^2-a1.^2) ;

            Sigma = (1-R1) .* (1-R2) .* (a2.*obj.P./(2.*pi.*k2)) .* (exp(-(del.*obj.PumpRadius).^2./8)./(B2.^2-a2.^2)) .* exp(-a1.*L1) ;

            A = -(  (1-g).*(b.*gamma.*(1-n.*v)-v).*exp(-a1.*L1-B1.*L1) + (g+v).*(1+b.*gamma.*(1+n))).*Gamma./H ...
                + ( (1-g).*b.*(vv-s+gamma.*(1+vv)).*exp(-B1.*L1)).*Sigma./H ;

            B = -( (1+g).*(b.*gamma.*(1-n.*v)-v).*exp(-a1.*L1-B1.*L1)+(g+v).*(1-b.*gamma.*(1-n)).*exp(-2.*B1.*L1)).*Gamma./H +...
                ( (1+g).*b.*(vv-s+gamma.*(1+vv)).*exp(-B1.*L1) ).*Sigma./H ;

            obj.EE = Gamma+A+B ;
            obj.B0 = B0 ;
            obj.Lth0 = Lth0 ;
        end
        function region0_GenerateMeshes(obj)
            % Generates a variety of meshes needed to calculate the
            % deflection
            % Define the range of x and y values to calculate Temperature/deflection
            obj.y = gpuArray(0:5*10^-6:2500*10^-6) ;
            obj.x = obj.xdata-obj.xshift ;

            % Generate ray sampling mesh. This defines the light rays where
            % deflection will be calculated. The end result is a sum of these rays.
            obj.xwidth = obj.ProbeRadius*4 ;
            obj.zwidth = obj.Lth0/4 ;  % Perform computation in steps of 1/8 thermal lengths
            if length(obj.x)>1
                obj.xstep = (obj.x(2)-obj.x(1))/3 ;
            else
                obj.xstep = 3e-6 ;
            end

            obj.x_raySample = (obj.x(1)-2*obj.xwidth):obj.xstep:(obj.x(end)+2*obj.xwidth) ;
            obj.z_raySample = -obj.zwidth-obj.zwidth*4*6 :obj.zwidth:-1e-8 ;

            % Precalculate dT0/dz for a range of r values, at several z heights
            rmax = sqrt( max(obj.y)^2 + max(obj.x)^2)*1.2 ;
            obj.r_precalc = gpuArray(0:rmax/1000:rmax) ;
            obj.T0dzOut_precalc = gpuArray(zeros(length(obj.r_precalc),length(obj.z_raySample))) ;
        end
        function Temp = region0_T0dz(obj, z_in)
            % This function computes the derivative of the temperature profile with
            % respect to the z direction ( normal to the surface). This is an important
            % quantity when calculating the delfection of a ray due to the gradient in
            % refractive index.
            
            delmatrix = obj.del'*ones(1,length(obj.r_precalc)) ;

            T0int =  (delmatrix).*besselj(0,obj.del'*obj.r_precalc).*(obj.EE'*ones(1,length(obj.r_precalc))).*(obj.B0'*ones(1,length(obj.r_precalc))).*exp(obj.B0'*ones(1,length(obj.r_precalc)).*z_in) ;

            % Do simpsons rule
            step = obj.del(3)-obj.del(1) ;
            index_vector = 1:((length(obj.del)-1)/2) ;
            Temp = sum((step/6).*(     T0int(2.*index_vector-1 ,:) ...
                       +4.*T0int(2.*index_vector ,:) ...
                       +   T0int(2.*index_vector+1 ,:)) ,1) ; 
                   
        end
        function region0_T0Calculate(obj)
            % This function determines the temperature profile for a number
            % of given z and x probe beam positions
            for jj=1:length(obj.z_raySample)
                obj.T0dzOut_precalc(:,jj) = obj.region0_T0dz(obj.z_raySample(jj)) ;
            end
        end
        function region0_CalcNonWeightedDeflection(obj)
            % Calculates the deflection profile without weighting by the
            % probe beam intensity.
            PsiNormal_gpu = gpuArray(zeros(length(obj.x_raySample),length(obj.z_raySample))) ;

            %%%%% This part =0.0444 sec %%%%%%%%

            r_matlab = sqrt(   (obj.x_raySample'*ones(1,length(obj.y))).^2 + ...
                               (ones(length(obj.x_raySample),1)*obj.y).^2      ) ;
            r_gpu = gpuArray(r_matlab) ;

            for zz=1:length(obj.z_raySample)

                % Interpolation is used to determine deflection values at
                % intermediate points
                T0dzOut = interp1(obj.r_precalc,obj.T0dzOut_precalc(:,zz),r_gpu) ;

                PsiNormal_gpu(:,zz) = 2.*obj.n.*obj.dNdT.*trapz(obj.y,T0dzOut,2) ;

            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            obj.PsiNormal = gather(PsiNormal_gpu)' ;
        end
        function region0_CalcFinalDeflection(obj)
            % Calculates the final, effective beam deflection by weighting
            % the deflection of each ray by the probe beam intensity
            
            
            % Reshape matrix and define Z,X meshes
            PsiNormalMesh = obj.PsiNormal ;
            [X,Z] = meshgrid(obj.x_raySample,obj.z_raySample) ;
            obj.PsiFinal = gpuArray(zeros(1,length(obj.x))) ;

            z_raySample_probe = obj.z_raySample ;
            xresolution = 40 ; % keep as an even number. This is the number of rays
                               % evaluated in the x direction for a single probe beam
                               % bundle
            
            % probe beam std to use for weigting calculation                   
            probeSigma = obj.ProbeRadius/( 2*sqrt(2*log(2))) ;
            for ii=1:length(obj.x)

            % This section will find the net deflection of a bundle of rays by
            % summing and weigting by the intensity profile of the probe beam

            x_raySample_probe = obj.x(ii)-obj.xwidth:2*obj.xwidth/(xresolution):obj.x(ii)+obj.xwidth ;

            [Xprobe,Zprobe] = meshgrid(x_raySample_probe,z_raySample_probe) ;
            
            Xlb = x_raySample_probe-obj.xwidth/(xresolution) ;
            Xub = x_raySample_probe+obj.xwidth/(xresolution) ;
            
            Zlb = z_raySample_probe-obj.zwidth/2 ;
            Zub = z_raySample_probe+obj.zwidth/2 ;  
            
            
            [Xlb_mesh, Zlb_mesh] = meshgrid(Xlb, Zlb) ;
            [Xub_mesh, Zub_mesh] = meshgrid(Xub, Zub) ;  
            
            obj.LB = [reshape(Xlb_mesh,[],1) reshape(Zlb_mesh,[],1)];  
            obj.UB = [reshape(Xub_mesh,[],1) reshape(Zub_mesh,[],1)];  

            PsiNormalProbe = interp2(X,Z,PsiNormalMesh,Xprobe,Zprobe) ;
            
            Sigma = [probeSigma.^2 probeSigma.^2] ;
            mu = [obj.x(ii) obj.z] ;
            obj.Weights = reshape(mvncdf(obj.LB,obj.UB,mu,Sigma),size(Xprobe))./mvncdf([-Inf -Inf],[Inf 0],mu,Sigma)  ;
           

            obj.PsiFinal(ii) =  sum(sum(PsiNormalProbe.* obj.Weights)) ;
            end

        end
        
        %Backside methods (experimental)
        
        function region3_CalculateCoeff(obj)
            % This function calculates the coefficents nessicary to compute the temeperature
            % profile of a bi-layer sample under periodic gaussian illumination.
            % 
            % More information can be found in this paper,
            % 
            % M. Reichling and H. Grönbeck, “Harmonic heat flow in isotropic layered 
            % systems and its use for thin film thermal conductivity measurements,” 
            % Journal of Applied Physics, vol. 75, no. 4, pp. 1914–1922, Feb. 1994.

            % Region 0, Fluid

            k0   = obj.Region0_k0 ;
            rho0 = obj.Region0_rho0 ;
            c0   = obj.Region0_c0 ;
            R1   = obj.Region0_R1 ;

            % Region 1, Thin Film

            k1   = obj.Region1_k1 ;
            rho1 = obj.Region1_rho1 ;
            c1   = obj.Region1_c1 ;
            a1   = obj.Region1_a1 ;
            L1   = obj.Region1_L1 ;
            R2   = obj.Region1_R2 ;

            % Region 2, Substrate

            k2   = obj.Region2_k2 ;
            rho2 = obj.Region2_rho2 ;
            c2   = obj.Region2_c2 ;
            a2   = obj.Region2_a2 ;
            L2   = obj.Region2_L2 ;
            Rth  = obj.Region2_Rth ;

            % Region 3, Fluid

            k3   = obj.Region3_k3 ;
            rho3 = obj.Region3_rho3 ;
            c3   = obj.Region3_c3 ;

            f = obj.f ; %#ok<*PROP>
            del = obj.del ;

            % Calculate
            Lth0 = sqrt( k0./(pi .* f .* rho0 .* c0)) ; % Checked
            Lth1 = sqrt( k1./(pi .* f .* rho1 .* c1)) ; % Checked
            Lth2 = sqrt( k2./(pi .* f .* rho2 .* c2)) ; % Checked
            Lth3 = sqrt( k3./(pi .* f .* rho3 .* c3)) ; % Checked

            B0 = sqrt( del.^2 + 2i./Lth0.^2 ) ; % Checked
            B1 = sqrt( del.^2 + 2i./Lth1.^2 ) ; % Checked
            B2 = sqrt( del.^2 + 2i./Lth2.^2 ) ; % Checked
            B3 = sqrt( del.^2 + 2i./Lth3.^2 ) ; % Checked

            g = k0.*B0./(k1.*B1) ; % Checked
            b = k2.*B2./(k1.*B1) ; % Checked
            v = a1./B1 ; % Checked
            n = k1.*B1.*Rth ; % Checked
            s = a2./B2 ; % Checked
            d = k3.*B3./(k2.*B2) ; % Checked
            mu = ((1-d)./(1+d)) .* exp(-2.*B2.*L2) ; % Checked
            vv = ((s-d)./(1+d)) .* exp(-(a2+B2).*L2) ; % Checked
            gamma = (1-mu)./(1+mu) ; % Checked



            H = (g+1).*(1+b.*gamma.*(1+n))-(1-b.*gamma.*(1-n)).*(1-g).*exp(-2.*B1.*L1) ; % Checked

            Gamma = (1-R1).*(a1.*obj.P).*exp(-(del.*obj.PumpRadius).^2./8)./(2.*pi.*k1)./(B1.^2-a1.^2) ; % Checked

            Sigma = (1-R1) .* (1-R2) .* (a2.*obj.P./(2.*pi.*k2)) .* (exp(-(del.*obj.PumpRadius).^2./8)./(B2.^2-a2.^2)) .* exp(-a1.*L1) ; % Checked

            A = -(  (1-g).*(b.*gamma.*(1-n.*v)-v).*exp(-a1.*L1-B1.*L1) + (g+v).*(1+b.*gamma.*(1+n))).*Gamma./H ...
                + ( (1-g).*b.*(vv-s+gamma.*(1+vv)).*exp(-B1.*L1)).*Sigma./H ; % Checked

            B = -( (1+g).*(b.*gamma.*(1-n.*v)-v).*exp(-a1.*L1-B1.*L1)+(g+v).*(1-b.*gamma.*(1-n)).*exp(-2.*B1.*L1)).*Gamma./H +...
                ( (1+g).*b.*(vv-s+gamma.*(1+vv)).*exp(-B1.*L1) ).*Sigma./H ; % Checked
            
            L = (1./(b.*(1-mu))) .*  ...
                (v.*Gamma.*exp(-a1.*L1)+A.*exp(-B1.*L1)-B.*exp(B1.*L1) + ...
                b.*(vv-s).*Sigma) ; % Checked
            
            G = vv.*Sigma + mu.*L ; % Checked
            
            M = Sigma.*exp(-a2.*L2)+L.*exp(-B2.*L2)+G.*exp(B2.*L2) ; % Checked

            obj.EE = Gamma+A+B ; % Checked
            obj.B0 = B0 ;
            obj.Lth0 = Lth0 ;
            
            obj.M = M ;
            obj.Lth3 = Lth3 ;
            obj.B3 = B3 ;
        end
        function region3_GenerateMeshes(obj)
            % Define the range of x and y values to calculate Temperature/deflection
            obj.y = gpuArray(0:5*10^-6:2500*10^-6) ;
            obj.x = obj.xdata-obj.xshift ;

            % Generate ray sampling mesh. This defines the light rays where
            % deflection will be calculated. The end result is a sum of these rays.
            obj.xwidth = obj.ProbeRadius*4 ;
            obj.zwidth = obj.Lth0/4 ;  % Perform computation over 4 thermal lengths

            if length(obj.x)>1
                xstep = (obj.x(2)-obj.x(1))/3 ;
            else
                xstep = 3e-6 ;
            end

            obj.x_raySample = (obj.x(1)-2*obj.xwidth):xstep:(obj.x(end)+2*obj.xwidth) ;
            obj.z_raySample = (obj.Region1_L1+obj.Region2_L2+obj.zwidth):obj.zwidth:(obj.Region1_L1+obj.Region2_L2+obj.zwidth + obj.zwidth*4*6) ;

            % Precalculate dT0/dz for a range of r values, at several z heights
            rmax = sqrt( max(obj.y)^2 + max(obj.x)^2)*1.2 ;
            obj.r_precalc = gpuArray(0:rmax/1000:rmax) ;
            obj.T3dzOut_precalc = gpuArray(zeros(length(obj.r_precalc),length(obj.z_raySample))) ;
        end
        function Temp = region3_T3dz(obj, z_in)
            % This function computes the derivative of the temperature profile with
            % respect to the z direction ( normal to the surface). This is an important
            % quantity when calculating the delfection of a ray due to the gradient in
            % refractive index.

            delmatrix = obj.del'*ones(1,length(obj.r_precalc)) ;
            
            T0int =  (delmatrix).*besselj(0,obj.del'*obj.r_precalc).*(obj.M'*ones(1,length(obj.r_precalc))).*(-obj.B3'*ones(1,length(obj.r_precalc))).*exp(-obj.B3'*ones(1,length(obj.r_precalc)).*(z_in-obj.Region1_L1-obj.Region2_L2)) ;

            % Do simpsons rule
            step = obj.del(3)-obj.del(1) ;
            index_vector = 1:((length(obj.del)-1)/2) ;

            Temp = sum((step/6).*(     T0int(2.*index_vector-1 ,:) ...
                       +4.*T0int(2.*index_vector ,:) ...
                       +   T0int(2.*index_vector+1 ,:)) ,1) ;
        end
        function region3_T3Calculate(obj)
            % This function determines the temperature profile
            for jj=1:length(obj.z_raySample)
                obj.T3dzOut_precalc(:,jj) = obj.region3_T3dz(obj.z_raySample(jj)) ;
            end
        end
        function region3_CalcNonWeightedDeflection(obj)
            
            PsiNormal_gpu = gpuArray(zeros(length(obj.x_raySample),length(obj.z_raySample))) ;

            %%%%% This part =0.0444 sec %%%%%%%%

            r_matlab = sqrt(   (obj.x_raySample'*ones(1,length(obj.y))).^2 + ...
                               (ones(length(obj.x_raySample),1)*obj.y).^2      ) ;
            r_gpu = gpuArray(r_matlab) ;

            for zz=1:length(obj.z_raySample)

                % Interpolation is used to determine deflection values at
                % intermediate points
                T3dzOut = interp1(obj.r_precalc,obj.T3dzOut_precalc(:,zz),r_gpu) ;

                % Numerical integration along y to determine deflection
                PsiNormal_gpu(:,zz) = 2.*obj.n.*obj.dNdT.*trapz(obj.y,T3dzOut,2) ;

            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            obj.PsiNormal = gather(PsiNormal_gpu)' ;
        end
        function region3_CalcFinalDeflection(obj)
            % Reshape matrix and define Z,X meshes
            PsiNormalMesh = obj.PsiNormal ;
            [X,Z] = meshgrid(obj.x_raySample,obj.z_raySample) ;
            obj.PsiFinal = gpuArray(zeros(1,length(obj.x))) ;

            z_raySample_probe = obj.z_raySample ;
            xresolution = 40 ; % keep as an even number. This is the number of rays
                               % evaluated in the x direction for a single probe beam
                               % bundle
            
            % probe beam std to use for weigting calculation                   
            probeSigma = obj.ProbeRadius/( 2*sqrt(2*log(2))) ;
            for ii=1:length(obj.x)

            % This section will find the net deflection of a bundle of rays by
            % summing and weigting by the intensity profile of the probe beam

            x_raySample_probe = obj.x(ii)-obj.xwidth:2*obj.xwidth/(xresolution):obj.x(ii)+obj.xwidth ;

            [Xprobe,Zprobe] = meshgrid(x_raySample_probe,z_raySample_probe) ;
            
            Xlb = x_raySample_probe-obj.xwidth/(xresolution) ;
            Xub = x_raySample_probe+obj.xwidth/(xresolution) ;
            
            Zlb = z_raySample_probe-obj.zwidth/2 ;
            Zub = z_raySample_probe+obj.zwidth/2 ;  
            
            
            [Xlb_mesh, Zlb_mesh] = meshgrid(Xlb, Zlb) ;
            [Xub_mesh, Zub_mesh] = meshgrid(Xub, Zub) ;  
            
            obj.LB = [reshape(Xlb_mesh,[],1) reshape(Zlb_mesh,[],1)];  
            obj.UB = [reshape(Xub_mesh,[],1) reshape(Zub_mesh,[],1)];  

            PsiNormalProbe = interp2(X,Z,PsiNormalMesh,Xprobe,Zprobe) ;
            
            Sigma = [probeSigma.^2 probeSigma.^2] ;
            mu = [obj.x(ii) obj.z] ;
            obj.Weights = reshape(mvncdf(obj.LB,obj.UB,mu,Sigma),size(Xprobe))./mvncdf([-Inf -Inf],[Inf 0],mu,Sigma)  ;
           

            obj.PsiFinal(ii) =  sum(sum(PsiNormalProbe.* obj.Weights)) ;
            end

        end


    end
    methods (Access = public) % Define setters/getters
        function Simulate(obj)
            % This is the primiary method used to run the calculation
            % the internal methods used in the calculation are optimzed
            % for calculation over the x position variable
            % Calculations over other parameters will nessicarilly be slow
            
            if obj.z < 0 
                % Frontside calc
                obj.region0_CalculateCoeff() % 0.004060
                obj.region0_GenerateMeshes() % 0.000790
                obj.region0_T0Calculate() % 1.3 sec
                obj.region0_CalcNonWeightedDeflection() % 0.048410 seconds
                obj.region0_CalcFinalDeflection() % 0.105022 seconds
            elseif obj.z > (obj.Region1_L1+obj.Region2_L2)
                % Backside calc
                obj.region3_CalculateCoeff()
                obj.region3_GenerateMeshes()
                obj.region3_T3Calculate()
                obj.region3_CalcNonWeightedDeflection
                obj.region3_CalcFinalDeflection()
            else
                error('Invalid z value')
            end
        end
        function out = get_norm_amplitude(obj)
            % Get results from simulation
            % Retrives normalized amplitude as a 1d array
            out(:,1) = abs(gather(obj.PsiFinal'.*exp(1i*obj.PhaseShift)))./max(abs(gather(obj.PsiFinal'.*exp(1i*obj.PhaseShift)))) ;
        end
        function out = get_amplitude(obj)
            % Get results from simulation
            % Retrives amplitude as a 1d array
            out(:,1) = abs(gather(obj.PsiFinal'.*exp(1i*obj.PhaseShift))) ;
        end
        function out = get_norm_phase(obj)
            % Output corrisponding normalized Phase image
            [~, I] = max(abs(gather(obj.PsiFinal'.*exp(1i*obj.PhaseShift)))) ;
            temp_data = unwrap(angle(gather(obj.PsiFinal'))) ;
            ref_angle = temp_data(I) ;
            out(:,1) = temp_data-ref_angle ;
        end
        function out = get_phase(obj)
            % Output corrisponding Phase
            out = angle(gather(obj.PsiFinal'.*exp(1i*obj.PhaseShift))) ;
        end
        function out = get_complex_amplitude(obj)
            % Get results from simulation
            % Retrives complex amplitude as a 1d array
            out(:,1) = gather(obj.PsiFinal'.*exp(1i*obj.PhaseShift)) ;
        end
        function set_del_regions(obj, regions_in)
            % Set number of regions to compute temperature itegral over
            % More will give more accurate results, but more computation
            % time.
            obj.n_regions = regions_in ;
            obj.n_points = 2.*obj.n_regions +1 ;
            obj.del = gpuArray( ((1:obj.n_points)-1).*(obj.del_max/(obj.n_points-1))+obj.del_min) ;            
        end
        function set_del_max(obj,max_in)
            obj.del_max = max_in ;
            obj.n_points = 2.*obj.n_regions +1 ;
            obj.del = gpuArray( ((1:obj.n_points)-1).*(obj.del_max/(obj.n_points-1))+obj.del_min) ;                
        end
    end
    methods % Define constructor
        function obj = PDSExper3
            % Constructs PDSExpr object
            % Initializes GPU array. Can take some time.
            obj.n_points = 2.*obj.n_regions +1 ;
            obj.del = gpuArray( ((1:obj.n_points)-1).*(obj.del_max/(obj.n_points-1))+obj.del_min) ;
        end
    end
    methods % Misc helper functions
    % Define quick plotting methods
        function plot_real(obj, data_in)
            % Method for quick plotting exepremental data
            % and store theoretical model
            % Only plots real components.
            h = figure
            map = [0, 0, 1
                   1, 0, 0
                   0, 1, 0
                   0, 1, 1
                   1, 0, 1
                   1, 1, 1];
            hold all
            for ii = (1:length(data_in))
            obj.xdata = data_in(ii).data(:,1) ;
            obj.f = data_in(ii).freq ;
            obj.Simulate()

            plot(data_in(ii).data(:,1),real(data_in(ii).data(:,4)),'color', map(ii,:),'Marker','o','LineStyle','none')
            
            plot(data_in(ii).data(:,1),real(obj.get_complex_amplitude()),'color', map(ii,:),'LineStyle','-','LineWidth',2)
            xlabel('Beam displacement (m)')
            ylabel('Beam deflection')
            end
        end
        function plot_imag(obj, data_in)
            % Method for quick plotting exepremental data
            % and store theoretical model
            % Only plots imaginary components.
            h = figure
            map = [0, 0, 1
                   1, 0, 0
                   0, 1, 0
                   0, 1, 1
                   1, 0, 1
                   1, 1, 1];
            hold all
            for ii = (1:length(data_in))
            obj.xdata = data_in(ii).data(:,1) ;
            obj.f = data_in(ii).freq ;
            obj.Simulate()

            plot(data_in(ii).data(:,1),imag(data_in(ii).data(:,4)),'color', map(ii,:),'Marker','o','LineStyle','none')
            
            plot(data_in(ii).data(:,1),imag(obj.get_complex_amplitude()),'color', map(ii,:),'LineStyle','-','LineWidth',2)
            xlabel('Beam displacement (m)')
            ylabel('Beam deflection')
            end
        end    
        function plot_amp(obj, data_in)
            % Method for quick plotting exepremental data
            % and store theoretical model
            % Plot amplitude
            map = [0, 0, 1
                   1, 0, 0
                   0, 1, 0
                   0, 1, 1
                   1, 0, 1
                   1, 1, 1];
            hold all
            for ii = (1:length(data_in))
            obj.xdata = data_in(ii).data(:,1) ;
            obj.f = data_in(ii).freq ;
            obj.Simulate()

            plot(data_in(ii).data(:,1),abs(data_in(ii).data(:,3)),'color', map(ii,:),'Marker','o','LineStyle','none')
            
            plot(data_in(ii).data(:,1),abs(obj.get_complex_amplitude()),'color', map(ii,:),'LineStyle','-','LineWidth',2)
            xlabel('Beam displacement (m)')
            ylabel('Beam deflection')
            end
        end  
        
        % Calculate residuals
        function Residual = resid_calc(obj, data_in, params)
            % Calculate squared residual given data/parameters
            % Used for fitting

            % Set parameters
            obj.Region1_k1 = params(1) ;
            obj.P = params(2) ;
            obj.z = params(3) ;
            obj.Region2_k2 = params(4) ;
            obj.Region2_Rth = params(5) ;
            obj.ProbeRadius = params(6) ;

            % Calculate residuals
            obj.xdata  = data_in(1).data(:,1) ;

            Residual=0 ;
            for ii=(1:length(data_in))
                obj.f=data_in(ii).freq ;

                % Note, the cost function is a weighted least squares term
                % The cost function is weighted by a gaussian
                % The mean and STD of the gassuan MUST be set properly such that
                % the relivant part of the data is fit to. IE the central peak
                % should be largely ignored. Instead the data should be mostly fit
                % to the decay

                obj.Simulate()
                Residual = Residual+ ...
                           sum( ((obj.get_amplitude() ...
                                - data_in(ii).data(:,3)  ).^2)) ; 
            end  

        end
        function Residual = logResid_calc(obj, data_in, params, cutoff)
            % Calculate the log squared residual
            % Used as a error function for fitting
            % Set parameters
            obj.Region1_k1 = params(1) ;
            obj.P = params(2) ;
            obj.z = params(3) ;
            obj.Region2_k2 = params(4) ;
            obj.Region2_Rth = params(5) ;

            % Calculate mask

            obj.xdata  = data_in(1).data(:,1) ;

            Residual=0 ;
            for ii=(1:length(data_in))
                obj.f=data_in(ii).freq ;
                mask = (data_in(ii).data(:,3) > cutoff) ;
                obj.xdata  = data_in(1).data(mask,1) ;

                % Note, the cost function is a weighted least squares term
                % The cost function is weighted by a gaussian
                % The mean and STD of the gassuan MUST be set properly such that
                % the relivant part of the data is fit to. IE the central peak
                % should be largely ignored. Instead the data should be mostly fit
                % to the decay

                obj.Simulate()
                Residual = Residual+ ...
                           sum( ((log(obj.get_amplitude()) ...
                                - log(data_in(ii).data(mask,3)) ).^2)) ; 
            end          
        end
        function y_final = mirror_correct(obj, data_in, ii, center)
            % Flips input data around center and averages
            % Used for correcting assymetry in measurements
            xdat = data_in(ii).data(:,1) ;
            ydat = data_in(ii).data(:,3) ;
            xdat2= -(xdat-center)+center ;
            ydat2 = interp1(xdat2, ydat, xdat) ;
            y_final = (ydat+ydat2)/2 ;
            y_final(isnan(y_final)) = 0 ;
        end

    
    end
end