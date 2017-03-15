classdef DataSet
    %% Define Properties
    properties (Constant)
    Nt = 570;
    Nd = 8;
    Ntr = 100;
    dt = 0.001;
    end
    
    properties
    Nn = 98
    Dir = {};    
    end
    
    methods
        %% Constructor
        function D = DataSet(data,N)  
            %Look for every direction
            for d = 1:D.Nd
                %Temp Struct Time
                temp = struct;
                S = zeros(D.Nn,D.Nt);
                X = zeros(2,D.Nt+1);
                V = zeros(2,D.Nt);  
                %Look for every trial
                for n = N(1):N(2)
                    %Local Data    
                    trial = data(n,d);
                    s = trial.spikes(:,1:D.Nt);
                    x = trial.handPos(1:2,1:D.Nt+1);
                    v = diff(x')';
                    %Spike, Position, Velocity,TrialMean
                    %Trial
                    S = S + s;
                    X = X + x;
                    V = V + v;
                end   
            %Average   
            Ntr = N(2) - N(1) + 1;   
            temp.Spikes = S/Ntr;
            temp.Position = X/Ntr;
            temp.Velocity = V/Ntr;
            %Append
            D.Dir{d} = temp;
            end
        end
        %% Operations
        %Convolution
        function D = Convolution(D,w,FieldIn,FieldOut)
            Nw = length(w);
            Nw_half = floor(Nw/2);
            for d = 1:D.Nd
                temp = [];
                    for i = 1:D.Nn
                        f = D.Dir{d}.(FieldIn);
                        %Convolution with Window
                        r = conv(w,f(i,:));
                        %Trim
                        temp(i,:) = r(Nw_half+1:end-Nw_half);
                    end
                D.Dir{d}.(FieldOut) = temp;
            end 
        end
        %% Edit
        %Eliminate Neurons
        function D = EliminateUnit(D,Index,Field)
            for d = 1:D.Nd
                D.Dir{d}.(Field) = removerows(D.Dir{d}.(Field),'ind',Index);
            end
            D.Nn = D.Nn - length(Index);
        end
        %Normalisation
        function D = BaseLineNormalisation(D,Field,Tcut)
            for d = 1:D.Nd
                F = D.Dir{d}.(Field);
                f = [];
                for i = 1:D.Nn
                    %Base Firing
                    B = mean(F(i,1:Tcut));
                    %Normalise
                    f(i,:) = (F(i,:)-B)/(max(F(i,:))- B); 
                end
                %Append
                D.Dir{d}.(Field) = f;
            end
        end
        %Displace
        function D = Lag(D,Field,Lag)
            for d = 1:D.Nd
                X = D.Dir{d}.(Field);
                [Nx,~] = size(X);
                Y = zeros(Nx,D.Nt);
                XLag = lagmatrix(X',-Lag)';
                Y(:,1:end-Lag) = XLag(:,1:end-Lag);
                D.Dir{d}.(Field) = Y;
            end
        end
        
        %% Output
        function [W,B,E] = GetPreDirection(D,a)
        %Preferred Direction
        temp = [];
            V = [];
            F = [];
            for d = 1:D.Nd
                v = D.Dir{d}.Velocity';
                V = [V;v];
                %Firing Matrix
                f = D.Dir{d}.FiringRate'; 
                F = [F;f];               
            end
            %Add Line Constant
            [Nv,~]  = size(V); 
            V = [ones(Nv,1),V];
            B = (V'*V+a(1)*ones(3,3))^-1*V'*F;
            %Error
            E = F-V*B;
            %Variance of the Residuals
            s = ones(1,D.Nn)./var(E);
            S = diag(s);
            W = pinv(B*S*B' + a(2)*ones(3,3))*B*S;
        end
        
        %% Test
        function [Ve] = MeanTest(W,Test)
            %Create Figure
            h1 = figure();
            hold on
            axis equal
            h2 = figure();
            subplot(2,4,1);
            %For Each Direction
            for d = 1:Test.Nd
                %Plot X and V
                F = Test.Dir{d}.FiringRate;
                X = zeros(3,Test.Nt+1);
                V = zeros(3,Test.Nt);
                    for t = 1:Test.Nt
                        %Regress
                        V(:,t) = W*(F(:,t));
                        X(:,t+1) = X(:,t) + V(:,t);    
                    end
                %Position
                figure(h1)
                hold on
                plot(X(2,:),X(3,:),'b')
                Xreal = Test.Dir{d}.Position;
                Xreal = Xreal - Xreal(:,1);
                plot(Xreal(1,:),Xreal(2,:),'r')
                %Velocity 
                Vreal = Test.Dir{d}.Velocity;
                Ve = Vreal - V(2:3,:);
                figure(h2)
                subplot(2,4,d)
                hold on
                l = plot(V');
                set(l(2),'Color','b')
                set(l(3),'Color','b','LineStyle','--')
                l = plot(Vreal');
                set(l(1),'Color','r')
                set(l(2),'Color','r','LineStyle','--')
            end 
        end
        
    end

    methods (Static)
        %Produce Unit Vector
        function V = UnitVector(V,Nn)
           for i = 1:Nn
               v = norm(V(:,i));
               V(:,i) = V(:,i)/v; 
           end
        end
        %Test
        function [Fe] = TrialTest(B,w,trial,N)
            Fe = [];
            h = waitbar(0,'Computing Noise Covariance');    
            for n = N(1):N(2)
                Test = DataSet(trial,[n n]);
                Test = EliminateUnit(Test,[38,49,52,76],'Spikes');
                Test = Convolution(Test,w,'Spikes','FiringRate');
                for d = 1:Test.Nd
                    V = [ones(Test.Nt,1), Test.Dir{d}.Velocity'];
                    F = Test.Dir{d}.FiringRate';
                    Fe = [Fe;F - V*B]; 
                end
                waitbar((n-N(1))/(N(2)-N(1)))
            end
            close(h);
        end   
    end
end