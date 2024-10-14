
%--------------------------------------------------------------------------
%Subject: numerical routine to solve flux flow in porous media
%Type of file: FUNCTION
%Programer: Fernando R. L. Contreras
%--------------------------------------------------------------------------
%Goals: this FUNCTION gets the value to boundary condition according the
%benchmark which intend to run.

%--------------------------------------------------------------------------
%Aditional comments:

%--------------------------------------------------------------------------

%Fill the matrix of permeability as a function of element to left and right
%of half edge evaluated. This function receives "kmap" and a feature of
%element which wants to know the permeability ("kfeature").
function [kmap] = PLUG_kfunction(kmap,h,MM)
%Define global parameters:
global elem centelem numcase coord;

%Choose the benchmark to attribute permeability.
switch numcase
    case 341
        [auxperm2,]=ferncodes_calcpermeab;
        for ii=1:size(centelem,1)
            kmap(ii,:)=[ii auxperm2(ii) 0 0 auxperm2(ii)];
        end
    case 342
        kmap(1,1:5) = [1 MM*kmap(1,2) MM*kmap(1,3) MM*kmap(1,4) MM*kmap(1,5)];
        elem(:,5)=1;
        %----------------------------------------------------------------------
        %Example 1.7: Axisymmetric Case (Anisotropic Medium). teta = pi/6
    case 333
        %Initialize "kaux"
        kaux = zeros(size(centelem,1),5);
        for ii = 1:size(centelem,1)

            kaux(ii,:) = [ii h(ii)*0.5 0 0 h(ii)*0.5];

        end  %End of FOR

        %Restore "kmap"
        kmap = kaux;
    case 331
        %Initialize "kaux"
        kaux = zeros(size(centelem,1),5);
        for ii = 1:size(centelem,1)

            kaux(ii,:) = [ii h(ii)*33.3 0 0 h(ii)*33.33];

        end  %End of FOR

        %Restore "kmap"
        kmap = kaux;
    case 1.7
        %Definition of "R" matrix
        %Initialization
        kk = [kmap(1,2) kmap(1,3); kmap(1,4) kmap(1,5)];
        %Definition of angle
        teta = pi/6;
        %Definition of ratation matrix
        R = [cos(teta) sin(teta); -sin(teta) cos(teta)];
        %Define the permeability to be used
        kk = R'*kk*R;
        %Build "kmap"
        kmap(1,:) = [1 kk(1,1) kk(1,2) kk(2,1) kk(2,2)];

        %----------------------------------------------------------------------
        %Example 1.8: Axisymmetric Case (Heterogeneous and Anisotropic Medium).
        %teta = pi/6
    case 1.8
        %Initialize a parameter
        epsilon = 1e-3;
        theta = pi/6;
        kk = [1000 0; 0 1];
        %Rotate the tensor in "theta"
        Krot = [cos(theta) -sin(theta); sin(theta) cos(theta)]*kk*...
            [cos(theta) sin(theta); -sin(theta) cos(theta)];

        %Initialize "kmap"
        kmap = zeros(size(centelem,1),5);
        for ii = 1:size(centelem,1)
            %Define "x" and "y"
            r = centelem(ii,1);
            z = centelem(ii,2);
            %Choose the tensor according "x" position
            if z > 0.5
                kmap (ii,:) = [ii Krot(1,1) 0 0 Krot(2,2)];
            else
                %Define "x1" and "y1"
                r1 = r + 1e-3;
                z1 = z + 1e-1;
                %Definition of permeability components
                kk(1,1) = ((z1^2) + epsilon*(r1^2));
                kk(1,2) = -(1 - epsilon)*(r1*z1);
                kk(2,1) = -(1 - epsilon)*(r1*z1);
                kk(2,2) = ((r1^2) + epsilon*(z1^2));
                %Build "kmap"
                kmap(ii,:) = [ii kk(1,1) 0 0 kk(2,2)];
            end  %End of IF
        end  %End of FOR

        %----------------------------------------------------------------------
        %Example 7.1: Adaptaded from FVCA 5 (Herbin and Hubert, 2008).
        %Case 3 - Oblique flow.
    case 7.1
        %Initialize "R":
        R = zeros(2);
        kk = [kmap(1,2) kmap(1,3); kmap(1,4) kmap(1,5)];
        %Fill "R"
        R(1,1) = cosd(40);
        R(1,2) = sind(40);
        R(2,1) = -R(1,2);
        R(2,2) = R(1,1);
        %Fill "k" turning the tensor
        kk = R'*kk*R;
        %Buld "kmap" again
        kmap = [1 kk(1,1) kk(1,2) kk(2,1) kk(2,2)];

        %----------------------------------------------------------------------
        %Example 7.2:  In this example there are two material with the first
        %isotropic and the second one orthotropic. Dirichlet's boundary
        %condition obtained from analitical solution (Drowniou and Le
        %Potier, 2011). Example 4.2.2 (Eq. 53)
    case 7.2
        %Initialize "kaux"
        kaux = zeros(size(centelem,1),5);
        for i = 1:size(centelem,1)
            %Define "x" and "y"
            x = centelem(i,1);
            if x <= 0.5
                kaux(i,:) = [i kmap(1,2:5)];
            else
                kaux(i,:) = [i kmap(2,2:5)];
            end  %End of IF
        end  %End of FOR

        %Restore "kmap"
        kmap = kaux;

        %----------------------------------------------------------------------
        %Example 15.1: Adaptaded from FVCA 5 (Herbin and Hubert, 2008).
        %It is the case 5 theirs (highly anisotropic), pp. 6
    case 15.1
        %Initialize a parameters
        delta = 1e-3;
        %Initialize "kmap"
        kmap = zeros(size(centelem,1),5);
        for i = 1:size(centelem,1)
            %Define "x" and "y"
            x = centelem(i,1);
            y = centelem(i,2);
            %Definition of permeability components
            k(1,1) = (delta*(x^2) + (y^2));
            k(1,2) = (delta - 1)*(x*y);
            k(2,1) = (delta - 1)*(x*y);
            k(2,2) = (delta*(y^2) + (x^2));

            k = (1/((x^2) + (y^2)))*k;

            %Build "kmap"
            kmap(i,:) = [i k(1,1) k(1,2) k(2,1) k(2,2)];
        end  %End of FOR

        %----------------------------------------------------------------------
        %Example 15.1: Adaptaded from FVCA 5 (Herbin and Hubert, 2008).
        %It is the case 5 theirs (highly anisotropic), pp. 6
    case 15.2
        %Initialize a parameters
        delta = 100;
        %Initialize "kmap"
        kmap = zeros(size(centelem,1),5);
        for i = 1:size(centelem,1)
            %Define "x" and "y"
            x = centelem(i,1);
            y = centelem(i,2);
            %Definition of permeability components
            k(1,1) = (delta*(x^2) + (y^2));
            k(1,2) = (delta - 1)*(x*y);
            k(2,1) = (delta - 1)*(x*y);
            k(2,2) = (delta*(y^2) + (x^2));

            k = (1/((x^2) + (y^2)))*k;

            %Build "kmap"
            kmap(i,:) = [i k(1,1) k(1,2) k(2,1) k(2,2)];
        end  %End of FOR

        %----------------------------------------------------------------------
        %Example 15.3: Adaptaded from FVCA 5 (Herbin and Hubert, 2008).
        %It is the case 5, pp. 6 (highly anisotropic) - MODIFIED by Le
        %Potier (Section 3.1, Eq. 21 e 22 from papaer "Finite volume scheme
        %satisfying maximum and minimum principles"). There is a sutil
        %difference between this example and the another two afore.
    case 15.3
        %Define parameters:
        epsilon = 1e-3;
        %Initialize "kmap"
        kmap = zeros(size(centelem,1),5);
        for i = 1:size(centelem,1)
            %Define "x1" and "y1"
            x1 = centelem(i,1) + epsilon;
            y1 = centelem(i,2) + epsilon;
            %Definition of permeability components
            k(1,1) = ((y1^2) + epsilon*(x1^2));
            k(1,2) = -(1 - epsilon)*(x1*y1);
            k(2,1) = -(1 - epsilon)*(x1*y1);
            k(2,2) = ((x1^2) + epsilon*(y1^2));
            %Build "kmap"
            kmap(i,:) = [i k(1,1) k(1,2) k(2,1) k(2,2)];
        end  %End of FOR

        %----------------------------------------------------------------------
        %Example 15.4: Adaptaded from FVCA 5 (Herbin and Hubert, 2008).
        %It is the case 5, pp. 6 (highly anisotropic) - MODIFIED by Le
        %Potier (Section 3.1, Eq. 21 e 22 from papaer "Finite volume scheme
        %satisfying maximum and minimum principles"). There is a sutil
        %difference between this example and the another two afore.
    case 15.4
        %Define parameters:
        epsilon = 1e-3;
        %Initialize "kmap"
        kmap = zeros(size(centelem,1),5);
        for i = 1:size(centelem,1)
            %Define "x1" and "y1"
            x1 = centelem(i,1) + epsilon;
            y1 = centelem(i,2) + epsilon;
            %Definition of permeability components
            k(1,1) = ((y1^2) + epsilon*(x1^2));
            k(1,2) = -(1 - epsilon)*(x1*y1);
            k(2,1) = -(1 - epsilon)*(x1*y1);
            k(2,2) = ((x1^2) + epsilon*(y1^2));
            %Build "kmap"
            kmap(i,:) = [i k(1,1) k(1,2) k(2,1) k(2,2)];
        end  %End of FOR

        %----------------------------------------------------------------------
        %Example 16:  In this example there are two material with the first
        %isotropic and the second one orthotropic. Dirichlet's boundary
        %condition obtained from analitical solution (Drowniou and Le
        %Potier, 2011). Example 4.2.1 (Eq. 51 and 52)
    case 16
        %Initialize "kaux"
        kaux = zeros(size(centelem,1),5);
        for i = 1:size(centelem,1)
            %Define "x" and "y"
            x = centelem(i,1);
            if x <= 0.5
                kaux(i,:) = [i kmap(1,2:5)];
            else
                kaux(i,:) = [i kmap(2,2:5)];
            end  %End of IF
        end  %End of FOR

        %Restore "kmap"
        kmap = kaux;

        %----------------------------------------------------------------------
        %Example 21: Lipnikov et al., 2007 (Example 1)
    case 21
        %Initialize a parameter
        epsilon = 5e-2;
        %Initialize "kmap"
        kmap = zeros(size(centelem,1),5);
        for i = 1:size(centelem,1)
            %Define "x" and "y"
            x = centelem(i,1);
            y = centelem(i,2);
            %Definition of permeability components
            k(1,1) = (y^2) + (epsilon*(x^2));
            k(1,2) = -(1 - epsilon)*x*y;
            k(2,1) = -(1 - epsilon)*x*y;
            k(2,2) = (epsilon*(y^2)) + (x^2);
            %Build "kmap"
            kmap(i,:) = [i k(1,1) k(1,2) k(2,1) k(2,2)];
        end  %End of FOR

        %----------------------------------------------------------------------
        %Example 21.1: Adapted from Lipnikov et al., 2007 (Example 1)
        %It changes the parameter "epsilon"
    case 21.1
        %Initialize a parameter
        epsilon = 5e-5;
        %Initialize "kmap"
        kmap = zeros(size(centelem,1),5);
        for i = 1:size(centelem,1)
            %Define "x" and "y"
            x = centelem(i,1);
            y = centelem(i,2);
            %Definition of permeability components
            k(1,1) = (y^2) + (epsilon*(x^2));
            k(1,2) = -(1 - epsilon)*x*y;
            k(2,1) = -(1 - epsilon)*x*y;
            k(2,2) = (epsilon*(y^2)) + (x^2);
            %Build "kmap"
            kmap(i,:) = [i k(1,1) k(1,2) k(2,1) k(2,2)];
        end  %End of FOR

        %----------------------------------------------------------------------
        %Example 22: Lipnikov et al., 2007 (Example 2)
    case 22
        %Initialize "R":
        R = zeros(2);
        %Define parameters:
        for i = 1:size(kmap,1)
            %Define parameters:
            material = kmap(i,1);
            k = [kmap(i,2) kmap(i,3); kmap(i,4) kmap(i,5)];
            %Choose the material number located in fifth column of "elem"
            switch material
                case 1
                    %Fill "R"
                    R(1,1) = cos(pi/6);
                    R(1,2) = sin(pi/6);
                    R(2,1) = -R(1,2);
                    R(2,2) = R(1,1);
                    %Fill "k" turning the tensor
                    k = R'*k*R;
                    %Build "kmap"
                    kmap(i,2:5) = [k(1,1) k(1,2) k(2,1) k(2,2)];
                case 2
                    %Fill "R"
                    R(1,1) = cos(-pi/6);
                    R(1,2) = sin(-pi/6);
                    R(2,1) = -R(1,2);
                    R(2,2) = R(1,1);
                    %Fill "k" turning the tensor
                    k = R'*k*R;
                    %Build "kmap"
                    kmap(i,2:5) = [k(1,1) k(1,2) k(2,1) k(2,2)];
                case 3
                    %Fill "R"
                    R(1,1) = cos(pi/6);
                    R(1,2) = sin(pi/6);
                    R(2,1) = -R(1,2);
                    R(2,2) = R(1,1);
                    %Fill "k" turning the tensor
                    k = R'*k*R;
                    %Build "kmap"
                    kmap(i,2:5) = [k(1,1) k(1,2) k(2,1) k(2,2)];
                case 4
                    %Fill "R"
                    R(1,1) = cos(-pi/6);
                    R(1,2) = sin(-pi/6);
                    R(2,1) = -R(1,2);
                    R(2,2) = R(1,1);
                    %Fill "k" turning the tensor
                    k = R'*k*R;
                    %Build "kmap"
                    kmap(i,2:5) = [k(1,1) k(1,2) k(2,1) k(2,2)];
            end  %End of SWITCH
        end  %End of FOR

        %----------------------------------------------------------------------
        %Example 23.1: Lipnikov et al., 2007 (Example 2), teta = pi/6
    case 23.1
        %Definition of "R" matrix
        %Initialization
        k = [kmap(1,2) kmap(1,3); kmap(1,4) kmap(1,5)];
        %Definition of angle
        teta = pi/6;
        %Definition of ratation matrix
        R = [cos(teta) sin(teta); -sin(teta) cos(teta)];
        %Define the permeability to be used
        k = inv(R)*k*R;
        %Build "kmap"
        kmap(1,:) = [1 k(1,1) k(1,2) k(2,1) k(2,2)];

        %----------------------------------------------------------------------
        %Example 23.2: Lipnikov et al., 2007 (Example 2), teta = pi/6
    case 23.2
        %Definition of "R" matrix
        %Initialization
        k = [kmap(1,2) kmap(1,3); kmap(1,4) kmap(1,5)];
        %Definition of angle
        teta = 5*pi/6;
        %Definition of ratation matrix
        R = [cos(teta) sin(teta); -sin(teta) cos(teta)];
        %Define the permeability to be used
        k = R'*k*R;
        %Build "kmap"
        kmap(1,:) = [1 k(1,1) k(1,2) k(2,1) k(2,2)];

        %----------------------------------------------------------------------
        %Example 23.3: Adapted from Lipnikov et al., 2007 (Example 1) and
        %Hubert and Herbin (2008), FVCA V (Example 5). It is a heterogeneous
        %and highly anisotropic media. The mesh is UNSTRUCTURED
    case 23.3
        %Initialize a parameter
        epsilon = 1e3;
        theta = 0.5*pi;
        k = [100 0; 0 0.01];
        %Rotate the tensor in "theta"
        Krot = [cos(theta) -sin(theta); sin(theta) cos(theta)]*k*...
            [cos(theta) sin(theta); -sin(theta) cos(theta)];

        %Initialize "kmap"
        kmap = zeros(size(centelem,1),5);
        for i = 1:size(centelem,1)
            %Define "x" and "y"
            x = centelem(i,1);
            y = centelem(i,2);
            %Choose the tensor according "x" position
            if x < 0.5
                kmap (i,:) = [i Krot(1,1) Krot(1,2) Krot(2,1) Krot(2,2)];
            else
                %Define "x1" and "y1"
                x1 = x + 1e-3;
                y1 = y + 1e-3;
                %Definition of permeability components
                k(1,1) = ((y1^2) + epsilon*(x1^2));
                k(1,2) = -(1 - epsilon)*(x1*y1);
                k(2,1) = -(1 - epsilon)*(x1*y1);
                k(2,2) = ((x1^2) + epsilon*(y1^2));
                %Build "kmap"
                kmap(i,:) = [i k(1,1) k(1,2) k(2,1) k(2,2)];
            end  %End of IF
        end  %End of FOR

        %----------------------------------------------------------------------
        %Example 23.4: Adapted from Lipnikov et al., 2007 (Example 1) and
        %Hubert and Herbin (2008), FVCA V (Example 5). It is a heterogeneous
        %and highly anisotropic media. The mesh is STRUCTURED
    case 23.4
        %Initialize a parameter
        epsilon = 1e-6;
        theta = 0.25*pi;
        k = [100 0; 0 0.01];
        %Rotate the tensor in "theta"
        Krot = [cos(theta) -sin(theta); sin(theta) cos(theta)]*k*...
            [cos(theta) sin(theta); -sin(theta) cos(theta)];

        %Initialize "kmap"
        kmap = zeros(size(centelem,1),5);
        for i = 1:size(centelem,1)
            %Define "x" and "y"
            x = centelem(i,1);
            y = centelem(i,2);
            %Choose the tensor according "x" position
            if x < 0.5
                kmap (i,:) = [i Krot(1,1) Krot(1,2) Krot(2,1) Krot(2,2)];
            else
                %Define "x1" and "y1"
                x1 = x + epsilon;
                y1 = y + epsilon;
                %Definition of permeability components
                k(1,1) = ((y1^2) + epsilon*(x1^2));
                k(1,2) = -(1 - epsilon)*(x1*y1);
                k(2,1) = -(1 - epsilon)*(x1*y1);
                k(2,2) = ((x1^2) + epsilon*(y1^2));
                %Build "kmap"
                kmap(i,:) = [i k(1,1) k(1,2) k(2,1) k(2,2)];
            end  %End of IF
        end  %End of FOR

        %----------------------------------------------------------------------
        %Example 23.5: Adapted from Lipnikov et al., 2007 (Example 1) and
        %Hubert and Herbin (2008), FVCA V (Example 5). It is a heterogeneous
        %and highly anisotropic media. The mesh is UNSTRUCTURED
    case 23.5
        %Initialize a parameter
        theta = 0.5*pi;
        k = [100 0; 0 0.01];
        %Rotate the tensor in "theta"
        Krot = [cos(theta) -sin(theta); sin(theta) cos(theta)]*k*...
            [cos(theta) sin(theta); -sin(theta) cos(theta)];

        %Initialize "kmap"
        kmap = zeros(size(centelem,1),5);
        for i = 1:size(centelem,1)
            %Define "x" and "y"
            x = centelem(i,1);
            %Choose the tensor according "x" position
            if x < 0.5
                kmap (i,:) = [i Krot(1,1) Krot(1,2) Krot(2,1) Krot(2,2)];
            else
                %Definition of "R" matrix
                %Initialization
                k = [1000 0; 0 1];
                %Definition of angle
                teta = pi/6;
                %Definition of ratation matrix
                R = [cos(teta) sin(teta); -sin(teta) cos(teta)];
                %Define the permeability to be used
                k = inv(R)*k*R;
                %Build "kmap"
                kmap(i,:) = [i k(1,1) k(1,2) k(2,1) k(2,2)];
            end  %End of IF
        end  %End of FOR

        %----------------------------------------------------------------------
        %Example 29: Adaptaded from Le Potier presentation
    case 29
        %Initialize a parameters
        e = 1e-3;
        %Initialize "kmap"
        kmap = zeros(size(centelem,1),5);
        for i = 1:size(centelem,1)
            %Define "x" and "y"
            x = centelem(i,1);
            y = centelem(i,2);
            %Define "x1" and "y1"
            x1 = x + e;
            y1 = y + e;
            %Definition of permeability components
            k(1,1) = ((y1^2) + e*(x1^2));
            k(1,2) = -(1 - e)*(x1*y1);
            k(2,1) = -(1 - e)*(x1*y1);
            k(2,2) = ((x1^2) + e*(y1^2));

            %Build "kmap"
            kmap(i,:) = [i k(1,1) k(1,2) k(2,1) k(2,2)];
        end  %End of FOR

        %----------------------------------------------------------------------
        %Example 31.2: Two-Phase Flow case. Adapted from Chueh et al., 2010
        %(Very High Permeability Distribution)
    case 31.2
        %         for i=1:size(centelem,1)
        %             k(i,:)=[i log(kmap(i,2)) 0 0 log(kmap(i,2)) ];
        %         end
        %        kmap=kmap;
        %Define number of the randomic values
        N = 40;
        %Define Randomic paramiter
        randcoord = getrandist;

        %Initialize "qsi" and "randcoord"
        qsi = zeros(N,1);

        %Initialize "kmap"
        kmap = zeros(size(centelem,1),5);

        %Swept all elements
        for i = 1:size(centelem,1)
            %Swept all randomic values in order define "randcoord"
            for irand = 1:N
                %Define "randcoord"
                qsi(irand) = exp(-(norm(centelem(i,1:2) - ...
                    randcoord(irand,:))/0.05)^2);
            end  %End of FOR

            %Definition of permeability constant
            kconst = min(max(sum(qsi),0.05),0.95);

            %Build "kmap"
            kmap(i,:) = [i kconst 0 0 kconst];
        end  %End of FOR

        %----------------------------------------------------------------------
        %Example 34.6: Two-Phase Flow case. Adapted from Chueh et al., 2010
        %(Very High Permeability Distribution for FIVE-SPOT Case)
    case 31.3
        kmap=kmap;

    case 34.6
        %        kthin = zeros((90)^2,5);
        %Define number of the randomic values
        N = 40;
        %Define Randomic paramiter
        randcoord = getrandist;
        %Initialize "qsi" and "randcoord"
        qsi = zeros(N,1);

        %Initialize "kmap"
        kmap = zeros(size(centelem,1),5);

        %         c = 1;
        %         m = 0;
        %         n = 90;
        %         o = 180;
        %Swept all elements
        for i = 1:size(centelem,1)
            %Swept all randomic values in order define "randcoord"
            for irand = 1:N
                %Define "randcoord"
                qsi(irand) = exp(-(norm(centelem(i,1:2) - ...
                    randcoord(irand,:))/0.05)^2);
            end  %End of FOR

            %Definition of permeability constant
            kconst = min(max(sum(qsi),0.05),0.95);

            %Build "kmap"
            kmap(i,:) = [i kconst 0 0 kconst];

            %             %Attribute to current column
            %             for j = 1:3
            %                 kthin(m + j,:) = [(m + j) kmap(i,2:5)];
            %                 kthin(n + j,:) = [(n + j) kmap(i,2:5)];
            %                 kthin(o + j,:) = [(o + j) kmap(i,2:5)];
            %             end
            %             if c == 30
            %                 m = m + 183;
            %                 n = n + 183;
            %                 o = o + 183;
            %                 c = 1;
            %             else
            %                 m = m + 3;
            %                 n = n + 3;
            %                 o = o + 3;
            %                 c = c + 1;
            %             end
        end  %End of FOR

        %         %Plot the file
        %         %Write table (Time (VPI), Oil Flow rate, Accumulated Oil and Water Cut)
        %         writeresult = fopen('C:\\Users\\Marcio\\Doutorado\\Programas\\kmap.dat','w');
        %         %Swept "results"
        %         for i = 1:size(kthin,1)
        %             %Write "result" according to "producelem" size
        %             %There is one producer well
        %              fprintf(writeresult,'%26.16E %26.16E %26.16E %26.16E %26.16E\r\n',...
        %                 kthin(i,1:5));
        %         end  %End of FOR (write data)
    case 34.7
        kmap=getchue;
    case 35
        for i=1:size(centelem,1)
            x=centelem(i,1);
            y=centelem(i,2);
            k0=1;
            s=1;
            m=3;
            l=3;
            r=k0*exp(sqrt(s)*cos(2*pi*m*x)*cos(2*pi*l*y));
            K(i,1:5) = [i r 0 0 r];
            %normKmap(i,1)= r;
            elem(i,5)=i;
            u=0;
        end
        kmap=K;
    case 36
        for ii=1:size(centelem,1)
            yy=centelem(ii,2);
            if yy<0.25
                K(ii,1:5) = [ii 1 0 0 1];
            elseif 0.25<yy && yy<0.5
                K(ii,1:5) = [ii 2 0 0 2];

            elseif 0.5<yy && yy<0.75
                K(ii,1:5) = [ii 5 0 0 5];
            else
                K(ii,1:5) = [ii 10 0 0 10];
            end
            u=0;
        end
        kmap=K;

    case 241
        %K(1,1:5) = [1 10 0 0 10]; % see Nilson's dissertation, when r=1
        %Initialization
        k = [10 0; 0 0.01];
        %Definition of angle
        %teta=0;
        teta = pi/7.2;
        %teta=pi/6;
        %Definition of ratation matrix
        R = [cos(teta) sin(teta); -sin(teta) cos(teta)];
        %Define the permeability to be used
        k = inv(R)*k*R;
        K(1,1:5) = [1 k(1,1) k(1,2) k(2,1) k(2,2)];
        %K(1,1:5) = [1 k(1,1) k(1,2) k(2,1) k(2,2)];% see Nilson's dissertation, when r=100
        kmap=K;
    case 245
        for i=1:size(centelem,1)
            x=centelem(i,1);
            y=centelem(i,2);
            k0=1;
            s=4;
            m=3;
            l=3;
            r=k0*exp(sqrt(s)*cos(2*pi*m*x)*cos(2*pi*l*y));
            K(i,1:5) = [i r 0 0 r];
            %normKmap(i,1)= r;
            elem(i,5)=i;
            u=0;
        end
        kmap=K;
    case 247
        for i=1:size(centelem,1)
            K(i,1:5) = [i kmap(i,1) 0 0 kmap(i,1)];
            %normKmap(i,1)= r;
            elem(i,5)=i;
        end
        kmap=K;
    case 249
        for i=1:size(centelem,1)
            K(i,1:5) = [i kmap(i,1) 0 0 kmap(i,1)];
            %normKmap(i,1)= r;
            elem(i,5)=i;
        end
        kmap=K;
    case 250
        for i=1:size(centelem,1)
            K(i,1:5) = [i kmap(i,1) 0 0 kmap(i,1)];
            %normKmap(i,1)= r;
            elem(i,5)=i;
        end
        kmap=K;
    case 248
        kmap=kmap;
    case 330
        % Permeability tensor for the case I, article 2023:
        % A Local Grid-Refined Numerical Groundwater Model
        % Based on the Vertex-centred Finite-Volume Method
        kmap(1,1:5) = MM*kmap;
        elem(:,5)=1;
    case 332
        % Permeability tensor for the case III, article 2023:
        % A Local Grid-Refined Numerical Groundwater Model
        % Based on the Vertex-centred Finite-Volume Method
        kmap(1,1:5) = [1 MM*kmap(1,2) MM*kmap(1,3) MM*kmap(1,4) MM*kmap(1,5)];
        elem(:,5)=1;
    case 334
        kmap(1,1:5) = [1 MM*kmap(1,2) MM*kmap(1,3) MM*kmap(1,4) MM*kmap(1,5)];
        elem(:,5)=1;
    case 335
        kmap(1,1:5) = [1 MM*kmap(1,2) MM*kmap(1,3) MM*kmap(1,4) MM*kmap(1,5)];
        elem(:,5)=1;
        %------------------------------------------------------------------
    case 336
        for j=1:size(centelem,1)
            if centelem(j,1)<0.5 ||centelem(j,1)==0.5
                K(j,1:5) = [j 1 0.5 0.5 1];
                elem(j,5)=j;
            else
                K(j,1:5) = [j 10 2 2 100];
                elem(j,5)=j;
            end
        end
        kmap=K;
    case 337
        kmap(1,1:5) = [1 MM*kmap(1,2) MM*kmap(1,3) MM*kmap(1,4) MM*kmap(1,5)];
        elem(:,5)=1;
    case 338
        kmap(1,1:5) = [1 MM*kmap(1,2) MM*kmap(1,3) MM*kmap(1,4) MM*kmap(1,5)];
        elem(:,5)=1;
    case 347
        %kmap(1,1:5) = [1 MM*kmap(1,2) MM*kmap(1,3) MM*kmap(1,4) MM*kmap(1,5)];
        kaux = zeros(size(centelem,1),5);
        for ii = 1:size(centelem,1)

            kaux(ii,:) = [ii h(ii)*9.26 0 0 h(ii)*9.26];
            elem(ii,5)=ii;
        end  %End of FOR
        kmap=kaux;
         
    case 43
        %Get the permeability field
        kmap = getnikitinperm(centelem);

        %----------------------------------------------------------------------
        %Example 43.2: Two-Phase Flow case. Adapted from Jizou and Riviere
        %8 litologies.
    case 43.2
        %Initialize "kmap"
        kmap = zeros(size(centelem,1),5);
        %Swept all elements:
        %         for i = 1:size(centelem,1)
        %             %Get the vertices:
        %             vertices = elem(i,1:4);
        %             %Get only the non zero values
        %             vertices = vertices(logical(vertices ~= 0));
        %             %Get the "y" coordinate of each vertex
        %             ycoordvtx = coord(vertices,2);
        %             %Calculate the maximum and minimum "y" coordinate
        %             ycoordmin = min(ycoordvtx);
        %             ycoordmax = max(ycoordvtx);
        %             %Define "x" and "y" (centroid of element)
        %             x = centelem(i,1);
        %             y = centelem(i,2);
        %             %Evaluate the position of "x" and "y"
        %             y_reg1 = -x + 0.25;
        %             y_reg2 = -x + 0.5;
        %             y_reg3 = -x + 0.75;
        %             y_reg4 = -x + 1;
        %             y_reg5 = -x + 1.25;
        %             y_reg6 = -x + 1.5;
        %             y_reg7 = -x + 1.75;
        %             %The element is in region 1 or 4
        %             if (x <= 0.25 && y <= y_reg1) || (x <= 0.25 && ...
        %                     y_reg1 >= ycoordmin && y_reg1 <= ycoordmax) ...
        %                     || (x > 0.75 && y >= y_reg7) || (x > 0.75 && ...
        %                     y_reg7 >= ycoordmin && y_reg7 <= ycoordmax)
        %                 %Definition of permeability components
        %                 k = [505 495; 495 505];
        %                 %The element is in region 2
        %             elseif (x <= 0.5 && y > y_reg1 && y < y_reg2) || ...
        %                     (y_reg2 >= ycoordmin && y_reg2 <= ycoordmax) || ...
        %                     (x > 0.5 && y < y_reg2)
        %                 %The element is in region 3
        %             elseif (x <= 0.5 && y > y_reg1 && y < y_reg2) || ...
        %                     (y_reg2 >= ycoordmin && y_reg2 <= ycoordmax) || ...
        %                     (x > 0.5 && y < y_reg2)
        %                 %Definition of permeability components
        %                 k = [1000 0; 0 10];
        %                 %The element is in region 3
        %             else
        %                 %Definition of permeability components
        %                 k = [10 0; 0 1000];
        %             end  %End of IF
        %             %Build "kmap"
        %             kmap(i,:) = [i k(1,1) k(1,2) k(2,1) k(2,2)];
        %         end  %End of FOR

        for i=1:size(centelem,1)
            if centelem(i,2)<=(0.251-centelem(i,1))
                elem(i,5)=i;

                k=[505 495; 495 505];

            elseif (0.251-centelem(i,1))<centelem(i,2) && centelem(i,2)<(0.51-centelem(i,1))
                elem(i,5)=i;

                k= [10 0; 0 1000];

            elseif (0.51-centelem(i,1))<=centelem(i,2) && centelem(i,2)<(0.751-centelem(i,1))
                elem(i,5)=i;

                k=[1000 0; 0 10];

            elseif (0.751-centelem(i,1))<=centelem(i,2) && centelem(i,2)<(0.99-centelem(i,1))
                elem(i,5)=i;

                k=[10 0; 0 1000];

            elseif (0.99-centelem(i,1))<=centelem(i,2) && centelem(i,2)<(1.249-centelem(i,1))
                elem(i,5)=i;

                k=[1000 0; 0 10];
            elseif (1.249-centelem(i,1))<=centelem(i,2) && centelem(i,2)<(1.49-centelem(i,1))
                elem(i,5)=i;

                k=[10 0; 0 1000];
            elseif (1.49-centelem(i,1))<=centelem(i,2) && centelem(i,2)<(1.749-centelem(i,1))
                elem(i,5)=i;

                k=[1000 0; 0 10];
            elseif (1.749-centelem(i,1))<=centelem(i,2)
                elem(i,5)=i;

                k=[505 495; 495 505];

            end
            %Build "kmap"
            kmap(i,:) = [i k(1,1) k(1,2) k(2,1) k(2,2)];
        end
        %----------------------------------------------------------------------
        %Example 44: Two-Phase Flow case. Adapted from Lipnikov et al., 2007
        %(Case 1)
    case 44
        %Initialize a parameter
        %         epsilon = 5e-3;
        %         %Initialize "kmap"
        %         kmap = zeros(size(centelem,1),5);
        %         for i = 1:size(centelem,1)
        %             %Define "x" and "y"
        %             x = centelem(i,1);
        %             y = centelem(i,2);
        %             %Definition of permeability components
        %             k(1,1) = (y^2) + (epsilon*(x^2));
        %             k(1,2) = -(1 - epsilon)*x*y;
        %             k(2,1) = -(1 - epsilon)*x*y;
        %             k(2,2) = (epsilon*(y^2)) + (x^2);
        %             %Build "kmap"
        %             kmap(i,:) = [i k(1,1) k(1,2) k(2,1) k(2,2)];
        %         end  %End of FOR
        k = [500 0; 0 1];
        %Fill "R"
        %theta=5*pi/6;
        %theta=67.5;
        theta=112;

        R(1,1) = cosd(theta);
        R(1,2) = sind(theta);
        R(2,1) = -R(1,2);
        R(2,2) = R(1,1);
        %Fill "k" turning the tensor
        A=inv(R);
        k = A*k*R;
        %Buld "kmap" again
        kmap = [1 k(1,1) k(1,2) k(2,1) k(2,2)];
        elem(:,5)=1;

        %----------------------------------------------------------------------
        %Example 46: Adapted from Lipnikov et al., 2007 (Example 1) and
        %Hubert and Herbin (2008), FVCA V (Example 5). It is a heterogeneous
        %and highly anisotropic media. The mesh is STRUCTURED.
        %OBS.: It is a two-phase vertion for case 23.4
    case 46
        %Initialize a parameter
        epsilon = 1e-3;
        theta = 0.25*pi;
        k = [1 0; 0 220];
        %         %Rotate the tensor in "theta"
        %         Krot = [cos(theta) -sin(theta); sin(theta) cos(theta)]*k*...
        %             [cos(theta) sin(theta); -sin(theta) cos(theta)];

        %Initialize "kmap"
        kmap = zeros(size(centelem,1),5);
        for i = 1:size(centelem,1)
            %Define "x" and "y"
            x = centelem(i,1);
            y = centelem(i,2);
            %Choose the tensor according "x" position
            %            if x < 0.5
            a11=x/(sqrt(x^2+y^2)); a12=y/(sqrt(x^2+y^2));
            %         %Rotate the tensor in "theta"
            Krot = [a11 -a12; a12 a11]*k*...
                [a11 a12; -a12 a11];

            kmap (i,:) = [i Krot(1,1) Krot(1,2) Krot(2,1) Krot(2,2)];
            %             else
            %                 %Define "x1" and "y1"
            %                 x1 = x + 1e-3;
            %                 y1 = y + 1e-3;
            %                 %Definition of permeability components
            %                 k(1,1) = ((y1^2) + epsilon*(x1^2));
            %                 k(1,2) = -(1 - epsilon)*(x1*y1);
            %                 k(2,1) = -(1 - epsilon)*(x1*y1);
            %                 k(2,2) = ((x1^2) + epsilon*(y1^2));
            %                 %Build "kmap"
            %                 kmap(i,:) = [i k(1,1) k(1,2) k(2,1) k(2,2)];
            %             end  %End of IF
        end  %End of FOR
    case 339
        for ii = 1:size(centelem,1)
            %Define "x" and "y"
            x = centelem(ii,1);
            if x<500 || x==500
                %  in this case considere T see pag. 23
                kmap (ii,:) = [ii 200 0 0 200];

            else
                kmap (ii,:) = [ii 20 0 0 20];
            end
        end
    case 380
        
        [riverresults]=xlsread('Teste_5.xlsx');
        kmap=zeros(size(riverresults,1),5);
        kmap(:,1)=1:size(riverresults,1);
        kmap(:,2)= riverresults(:,3);
        kmap(:,3)=0*riverresults(:,3);
        kmap(:,4)=0*riverresults(:,3);
        kmap(:,5)=riverresults(:,3);
        
        %[kmap]=conduchidraulica;
    case 380.1
        %kmap=kmap;
         [auxperm2,]=ferncodes_calcpermeab;
        for ii=1:size(centelem,1)
            kmap(ii,:)=[ii auxperm2(ii) 0 0 auxperm2(ii)];
        end
%         [riverresults]=xlsread('Teste_6.xlsx');
%         kmap=zeros(size(riverresults,1),5);
%         kmap(:,1)=1:size(riverresults,1);
%         kmap(:,2)= riverresults(:,3);
%         kmap(:,3)=0*riverresults(:,3);
%         kmap(:,4)=0*riverresults(:,3);
%         kmap(:,5)=riverresults(:,3);
%         elem(:,5)=1:size(elem,1);
end  %End of Switch

%--------------------------------------------------------------------------
%SPECIAL CASE: Kozdon et al. (2011)
%--------------------------------------------------------------------------
%Example 45: Two-Phase Flow case. Adapted from Kozdon et al., 2011.
%"Multidimensional upstream weighting for multiphase transport in
%porous media". Section 5.2
if numcase > 45 && numcase < 46
    %Initialize "kaux"
    kaux = zeros(size(centelem,1),5);
    %Define reference "r"
    r = (1 - (1/51))/2;
    for i = 1:size(centelem,1)
        %Define "rmesh"
        rmesh = norm(centelem(i,:));
        %Definition of permeability components
        %Internal permeability
        if rmesh < r
            kaux(i,:) = [i kmap(1,2:5)];
            %Boundary permeability
        elseif rmesh >= r
            kaux(i,:) = [i kmap(2,2:5)];
        end  %End of IF
    end  %End of FOR

    %Update "kmap"
    kmap = kaux;
    %Update "elem"
    elem(:,5) = 1:size(elem,1);
end  %End of IF (Kozdon)

%--------------------------------------------------------------------------
%Function "getnikitinperm"
%--------------------------------------------------------------------------
    function [kmap] = getnikitinperm(centelem)
        %Define global paramiters


        %Initialize "kmap"
        kmap = zeros(size(centelem,1),5);
        %Swept all elements:
        for i = 1:size(centelem,1)
            %Get the vertices:
            vertices = elem(i,1:4);
            %Get only the non zero values
            vertices = vertices(logical(vertices ~= 0));
            %Get the "y" coordinate of each vertex
            ycoordvtx = coord(vertices,2);
            %Calculate the maximum and minimum "y" coordinate
            ycoordmin = min(ycoordvtx);
            ycoordmax = max(ycoordvtx);
            %Define "x" and "y" (centroid of element)
            x = centelem(i,1);
            y = centelem(i,2);
            %Evaluate the position of "x" and "y"
            y_reg1 = -x + 0.5;  %-x + 0.75
            y_reg2 = -x + 1;
            y_reg3 = -x + 1.5;  %-x + 1.25
            %The element is in region 1 or 4
            if (x <= 0.5 && y <= y_reg1) || (x <= 0.5 && ...
                    y_reg1 >= ycoordmin && y_reg1 <= ycoordmax) ...
                    || (x > 0.5 && y >= y_reg3) || (x > 0.5 && ...
                    y_reg3 >= ycoordmin && y_reg3 <= ycoordmax)
                %Definition of permeability components
                k = [505 495; 495 505];
                %The element is in region 2
            elseif (x <= 0.5 && y > y_reg1 && y < y_reg2) || ...
                    (y_reg2 >= ycoordmin && y_reg2 <= ycoordmax) || ...
                    (x > 0.5 && y < y_reg2)
                %Definition of permeability components
                k = [1000 0; 0 10];
                %The element is in region 3
            else
                %Definition of permeability components
                k = [10 0; 0 1000];
            end  %End of IF
            %Build "kmap"
            kmap(i,:) = [i k(1,1) k(1,2) k(2,1) k(2,2)];
        end  %End of FOR
    end
end
