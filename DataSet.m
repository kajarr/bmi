classdef DataSet
    %% Define Properties
    properties (Constant)
    Nd = 8;
    Ntr = 50;
    NTrain = 1;
    Np = 10;
    end
    
    properties
    Nt = 700;
    Nn = 98
    NTest = 50;
    Train = {};
    Test = {};
    end
    
    %% Methods
    methods
        %% Constructor
        function [D] = DataSet(data,NSplit)  
            %% Training
            for d = 1:D.Nd
                %Temp Struct Time
                temp = struct;
                %Empty Array-Zero Pad
                S = zeros(D.Nn,D.Nt);
                %A = zeros(2,D.Nt-1);
                V = zeros(2,D.Nt);
                X = zeros(2,D.Nt+1);
                %Look for every trial
                for n = 1:NSplit
                    %Local Data    
                    trial = data(n,d);
                    %Spikes
                    s = trial.spikes(:,1:end);
                    [~,Nsample] = size(s);
                    %If to large
                    if Nsample > D.Nt
                        clear s
                        Nsample = D.Nt;
                        s = trial.spikes(:,1:Nsample);
                    end
                    S(:,1:Nsample) = S(:,1:Nsample) + s;
                    %Position-Velocity
                    x = trial.handPos(1:2,1:Nsample);
                    pad = [x(1,end)*ones(1,D.Nt+1-Nsample);x(2,end)*ones(1,D.Nt+1-Nsample)];
                    x = [x, pad];
                    X = X + x;
                    v =  diff(x')';
                    V = V + v;
                    clear s x v
                end    
            %Center
            Xo = zeros(2,D.Nt+1);
            Xo(1,:) = X(1,1);
            Xo(2,:) = X(2,1);
            %Average
            N = NSplit;
            temp.Spikes = S/N;
            temp.Position = (X-Xo)/N;
            temp.Velocity = V/N;
            %Polar
            [Xp,Vp] = DataSet.Polar(X);
            temp.PolarVelocity = Xp; 
            temp.PolarPosition = Vp;
            %Append
            D.Train{1,d} = temp;
            end
            
            %% Testing 
            for d = 1:D.Nd
                for n = NSplit+1:D.Ntr
                    %Temp Structure
                    temp = struct;
                    %Empty Array-Zero Pad
                    S = zeros(D.Nn,D.Nt);
                    V = zeros(2,D.Nt);
                    X = zeros(2,D.Nt+1);
                    %Local Data    
                    trial = data(n,d);
                    %Size
                    s = trial.spikes(:,1:end);
                    [~,Nsample] = size(s);
                    if Nsample > D.Nt
                        clear s
                        Nsample = D.Nt;
                        s = trial.spikes(:,1:Nsample);
                    end
                    %Padding
                    x = trial.handPos(1:2,1:Nsample);
                    pad = [x(1,end)*ones(1,D.Nt+1-Nsample);x(2,end)*ones(1,D.Nt+1-Nsample)];
                    x = [x, pad];
                    %Spikes
                    S(:,1:Nsample) = s;
                    %Position-Velocity
                    X(:,1:D.Nt+1) = x;
                    v =  diff(x')';
                    V = v;
                    %Center
                    Xo = zeros(2,D.Nt+1);
                    Xo(1,:) = X(1,1);
                    Xo(2,:) = X(2,1);
                    %Sort
                    temp.Spikes = S;
                    temp.Position = X-Xo;
                    temp.Velocity = V;
                    %Polar
                    [Xp,Vp] = DataSet.Polar(X);
                    temp.PolarVelocity = Xp; 
                    temp.PolarPosition = Vp;
                    %Append
                    D.Test{n - NSplit,d} = temp;
                    D.NTest = D.Ntr - NSplit;
                    clear s x v
                end    
            end
        end   
        
        %% Preprocessing
        %Convolution
        function D = Convolution(D,w)
            Nw = length(w);
            Nw_half = floor(Nw/2);
            Field = {'Train','Test','NTrain','NTest'};
            %Redirect
            for f = 1:2
                Set = Field{f};
                NSet = Field{f+2};
                for d = 1:D.Nd
                    for n = 1:D.(NSet)
                        temp = [];
                        for i = 1:D.Nn
                            %Get Data
                            s = D.(Set){n,d}.Spikes;
                            %Convolution 
                            r = conv(w,s(i,:));
                            %Trim
                            temp(i,:) = r(Nw_half+1:end-Nw_half+1);
                        end
                        D.(Set){n,d}.FiringRate = temp;
                    end
                end
            end
        end
        %Eliminate Neurons
        function D = EliminateUnit(D,Index,Field)
            for d = 1:D.Nd
                D.Train{1,d}.(Field) = removerows(D.Train{d}.(Field),'ind',Index);
                for n = 1:D.NTest
                D.Test{n,d}.(Field) = removerows(D.Test{n,d}.(Field),'ind',Index);
                end    
            end
            D.Nn = D.Nn - length(Index);
        end
        %Mean TimeStep
        function D = MeanTimeStep(D,dt)
            for d = 1:D.Nd
                x = D.Train{1,d}.Position;
                f = D.Train{1,d}.FiringRate;
                X = [];
                F = [];
                for t = dt:dt:D.Nt-dt
                    X = [X,x(:,t)];
                    F = [F,mean(f(:,t-dt/2:t+dt/2),2)];
                end
                X = [X,X(:,end)];
                V = diff(X')';
                [~,Nv] = size(V);
                Vn = [];
                for k = 1:Nv
                    Vn = [Vn,norm(V(:,k))];
                end
                D.Train{1,d}.FiringRate = F;
                D.Train{1,d}.Position = X;
                D.Train{1,d}.Velocity = V;
                D.Train{1,d}.VelocityNorm = Vn;
            end
            for n = 1:D.NTest
                for d = 1:D.Nd
                x = D.Test{n,d}.Position;
                f = D.Test{n,d}.FiringRate;
                X = [];
                F = [];
                for t = dt:dt:D.Nt-dt
                    X = [X,x(:,t)];
                    F = [F,mean(f(:,t-dt/2:t+dt/2),2)];
                end
                X = [X,X(:,end)];
                V = diff(X')';
                [~,Nv] = size(V);
                Vn = [];
                for k = 1:Nv
                    Vn = [Vn,norm(V(:,k))];
                end
                D.Test{n,d}.FiringRate = F; 
                D.Test{n,d}.Position = X;
                D.Test{n,d}.Velocity = V;
                D.Test{n,d}.VelocityNorm = Vn;   
                end
            end
            D.Nt = length(dt:dt:D.Nt-dt);
        end
        %PCA
        function [D,p] = PCA(D)
            %Collecting Firing Rate
            F = [];
            for d = 1:D.Nd
                f = D.Train{1,d}.FiringRate;
                F = [F,f];
            end
            %Covariance Plot
            Sf = cov(F');
            [P,L] = eig(Sf);
            L = diag(L);
            %Select Last Vector
            p = P(:,end-D.Np+1:end)';
            %Figures
            figure
            %Covariance Matrix
            imagesc(Sf);
            %Eigenvalue
            figure
            plot(L);
            %Create SubSpace
            figure
            for d = 1:D.Nd
                f = D.Train{1,d}.FiringRate;
                %Subspace 
                s = p*f;
                %Subspace p
                D.Train{1,d}.Subspace = s;
                subplot(2,4,d)
                plot(s')
                for n  = 1:D.NTest
                    f = D.Test{n,d}.FiringRate;
                    s = p*f;
                    D.Test{n,d}.Subspace = s;
                end
            end        
        end
        
        %% Regression
        %Preferred Direction
        function [B] = GetPreDirection(D,L)
            %Preferred Direction
            temp = [];
            V = [];
            F = [];
            for d = 1:D.Nd
                %Start at t = 320
                v = D.Train{1,d}.Velocity(:,14:end);
                vn = D.Train{1,d}.VelocityNorm(:,14:end);
                V = [V,[v;vn]];
                %Firing Matrix
                f = D.Train{1,d}.Subspace(:,14:end);
                %Select
                F = [F,f];               
            end
            %Add Line Constant
            [~,Nv]  = size(V); 
            V = [ones(1,Nv);V];
            B = F*V'*(V*V' + L*eye(4))^-1;
            figure
            hold on 
            for i = 1:D.Np
                plot3([0,B(i,2)],[0,B(i,3)],[0,B(i,4)])
            end
            figure
            e = mean(abs(F-B*V),2);
            bar(e);
        end
            
        %Auto-Regressive Model
        function [Aex] = AutoRegression(D,L)
            %State-Space
            Xold = [];
            Xnew = [];
            %Collecting
            P = [];
            V = [];
            for d = 1:D.Nd
                %Position-Velocity
                p = D.Train{1,d}.Position(:,1:end-1);
                P = [P,p];
                v = D.Train{1,d}.VelocityNorm;
                V = [V,v];
            end
            [X,x] = DataSet.Shift([P;V],2);
            %Regression
            A = X*x'*(x*x' + L*eye(6))^-1;
            %Zeros
            A(4:end,1:3) = eye(3);
            A(4:end,4:end) = zeros(3,3);
            %Extend
            Aex = zeros(7,7);
            Aex(1,1) = 1;
            Aex(2:end,2:end) = A;          
        end
        %Residuals
        function [Q,R] = Residual(D,A,B)
            EB = [];
            EA = [];
            h = waitbar(0,'Computing Residuals');
            for n = 1:D.NTrain
                for d = 1:D.Nd
                    %PD
                    V = [ones(1,D.Nt); D.Test{n,d}.Velocity;D.Test{n,d}.VelocityNorm];
                    F = D.Test{n,d}.Subspace;
                    EB = [EB,F - B*V]; 
                    p = [D.Train{1,d}.Position(:,1:end-1);D.Test{n,d}.VelocityNorm];
                    [Na,~] = size(A);
                    [X,x] = DataSet.Shift(p,(Na-1)/3);
                    X = [ones(1,D.Nt);X];
                    x = [ones(1,D.Nt);x];
                    EA = [EA,X - A*x];
                end
                waitbar(n/5);
            end
            close(h);
            %Covariance Matrix
            R = cov(EB');
            Q = cov(EA');
            figure
            subplot(1,2,1)
            imagesc(Q);
            subplot(1,2,2)
            imagesc(R);
        end
        %Linear Transformation
        function [Ac,bc] = LinearTransformation(D,model,L)
            %Model
            K = model.KFilter;
            zo = model.BaseFiring;
            % Reading
            R = [];
            P = [];
            for n = 1:D.NTrain
                for d = 1:D.Nd
                    p = D.Test{n,d}.Position(:,10:end-1); 
                    sub = D.Test{n,d}.Subspace(:,10:end);
                    xo = p(:,1);
                    X = [xo;xo];
                    r = [];
                    for t = 1:21
                        z = sub(:,t) - zo;
                        %KalmanFiltering
                        [K,X] = KalmanUpdate(K,X,z);
                        r = [r,[X(1);X(2)]];
                    end
                    hold on
%                     plot(r(1,:),r(2,:),'-ro')
%                     plot(p(1,:),p(2,:),'-bo')
                    R = [R,r];
                    P = [P,p];
                end
            end
            [~,Nr] = size(R);
            R1 = [ones(1,Nr);R];
            C = P*R1'*(R1*R1' + L*eye(3))^-1;
            bc = C(:,1);
            Ac = C(:,2:3);
%             for t = 1:Nr
%             p_est(:,t) = Ac*R(:,t) + bc;
%             end
%             plot(p_est(1,:),p_est(2,:),'-yo');
        end
        %% Plot
        %Tuning
        function [] = PlotTuning(D)
            axis equal
            hold on
            %Position
            for d = 1:D.Nd
                X = D.Train{1,d}.Position;
                plot(X(1,:),X(2,:));
            end
            figure
            %Velocity
            for d = 1:D.Nd
                subplot(2,4,d)
                plot(D.Train{1,d}.Velocity')
            end
%             figure
%             %Polar Velocity
%             for d = 1:D.Nd
%                 subplot(2,4,d)
%                 plot(D.Train{1,d}.PolarVelocity')
%             end   
            figure
            %Spike Train
            for d = 1:D.Nd
                subplot(2,4,d)
                imagesc(D.Train{1,d}.FiringRate)
            end
        end
    end
    
    %% Other Method
    methods (Static)
        %% Shift for Auto-Regression
        function [X,x] = Shift(y,n)
            Y  = {};
            for z = 1:n+1
            Y{z} = y(:,z:end-(n+1-z));
            end
            x = [];
            for z = 1:n
            x = [Y{z};x];
            end
            X = [];
            for z = 2:n+1
            X = [Y{z};X];
            end
        end
        %% Polar Coordinate Transform 
        function [Xp,Vp] = Polar(X)
           [phi,r] = cart2pol(X(1,:),X(2,:));
           %Continous Angle
           phi = unwrap(phi);
           %Velocity
           rdot = diff(r);
           w = diff(phi);
           %Velocity
           Xp = [r,phi];
           Vp = [rdot;r(2:end).*w];    
        end

        %% Windows
        %Gaussian kernel
        function f = gaussk(x, mu, s)
            ex = (-1/2)*((x - mu)/s).^2;
            ar = (s * sqrt(2*pi));
            f = (1/ar)*exp(ex);
        end
        %Exponential kernel
        function f = expk(x, mu, s)
            ex = (-sqrt(2))*abs((x-mu)/s);
            ar = (s * sqrt(2));
            f = (1/ar)*exp(ex);
        end
        %Causal kernel
        function hw = cauk(x,sp,a)
            q = ((a^2)*(x-sp));
            w = exp((-a).*(x-sp));
            for i = 1:length(q)
                f(i)=q(i)*w(i);
            end
            hw = zeros(size(f));
            for i = 1:length(f)
                if f(i) > 0.0
                    hw(i) = f(i);
                else
                    hw(i) = 0.0;
                end
            end
        end
    end
end
