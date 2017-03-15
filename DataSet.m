classdef DataSet
    %% Define Properties
    properties (Constant)
    Nt = 1000;
    Nd = 8;
    Ntr = 100;
    NTrain = 1;
    end
    
    properties
    Nn = 98
    NTest = 50;
    Train = {};
    Test = {};
    end
    
    %% Methods
    methods
        %% Constructor
        function D = DataSet(data,NSplit)  
            %% Training
            for d = 1:D.Nd
                %Temp Struct Time
                temp = struct;
                %Empty Array-Zero Pad
                S = zeros(D.Nn,D.Nt);
                A = zeros(2,D.Nt-1);
                V = zeros(2,D.Nt);
                X = zeros(2,D.Nt+1);
                %Look for every trial
                for n = 1:NSplit
                    %Local Data    
                    trial = data(n,d);
                    %Spikes
                    s = trial.spikes(:,1:end);
                    [~,Ns] = size(s);
                    S(:,1:Ns) = S(:,1:Ns) + s;
                    %Position-Velocity
                    x = trial.handPos(1:2,1:end);
                    [~,Nx] = size(x);
                    x = [x, x(:,end).*ones(2,D.Nt+1-Nx)];
                    X = X + x;
                    v =  diff(x')';
                    V = V + v;
                    a = diff(v')';
                    A = A + a;
                end    
            %Average
            N = NSplit;   
            temp.Spikes = S/N;
            temp.Position = X/N;
            temp.Velocity = V/N;
            temp.Acceleration = A/N;
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
                    A = zeros(2,D.Nt-1);
                    V = zeros(2,D.Nt);
                    X = zeros(2,D.Nt+1);
                    %Local Data    
                    trial = data(n,d);
                    %Size
                    s = trial.spikes(:,1:end);
                    x = trial.handPos(1:2,1:end);
                    [~,Ns] = size(s);
                    [~,Nx] = size(x);
                    %Spikes
                    S(:,1:Ns) = s;
                    %Position-Velocity
                    X(:,1:Nx) = x;
                    v =  diff(x')';
                    V(:,1:Nx-1) = v;
                    a = diff(v')';
                    A(:,1:Nx-2) = a;
                    %Append
                    temp.Spikes = S;
                    temp.Position = X;
                    temp.Velocity = V;
                    temp.Acceleration = A;
                    D.Test{n - NSplit,d} = temp;
                    D.NTest = D.Ntr - NSplit;
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
                            temp(i,:) = r(Nw_half+1:end-Nw_half);
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
        
        %% Regression
        %Preferred Direction
        function [B] = GetPreDirection(D,L)
        %Preferred Direction
            temp = [];
            V = [];
            F = [];
            for d = 1:D.Nd
                v = D.Train{1,d}.Velocity';
                V = [V;v];
                %Firing Matrix
                f = D.Train{1,d}.FiringRate'; 
                F = [F;f];               
            end
            %Add Line Constant
            [Nv,~]  = size(V); 
            V = [ones(Nv,1),V];
            B = (V'*V+L*ones(3,3))^-1*V'*F;
            figure
            hold on 
            for i = 1:D.Nn
                plot([0,B(2,i)],[0,B(3,i)])
            end
        end
        %Auto-Regressive Model
        function [A] = AutoRegression(D,L)
            %State-Space
            Xold = [];
            Xnew = [];
            %Collecting
            for d = 1:D.Nd
                %Position-Velocity
                bd = [290,650];
                p = D.Train{1,d}.Position(:,bd(1):bd(2));
                v = D.Train{1,d}.Velocity(:,bd(1):bd(2));
                %Concate and Shift
                %xold = [p(:,1:end-1);v(:,1:end-1)];
                %xnew = [p(:,2:end);v(:,2:end)]; 
                xold = v(:,1:end-1);
                xnew = v(:,2:end);
                Xold = [Xold,xold];
                Xnew = [Xnew,xnew];
            end
            %Regression
            A = Xnew*Xold'*(Xold*Xold' + L*eye(2))^-1;
        end
        %Residuals
        function [Q,R] = Residual(D,A,B)
            EB = [];
            EA = [];
            h = waitbar(0,'Computing Residuals');
            for n = 1:D.NTest
                for d = 1:D.Nd
                    %PD
                    V = [ones(D.Nt,1), D.Test{n,d}.Velocity'];
                    F = D.Test{n,d}.FiringRate';
                    EB = [EB;F - V*B]; 
                    %Velocity
                    v = D.Test{n,d}.Velocity;
                    x = D.Test{n,d}.Position(:,1:end-1);
                    %xold = [x(:,1:end-1);v(:,1:end-1)];
                    %xnew = [x(:,2:end);v(:,2:end)];
                    xold = [v(:,1:end-1)];
                    xnew = [v(:,2:end)]; 
                    EA = [EA,xnew-A*xold];
                end
                waitbar(n/D.NTest);
            end
            close(h);
            %Covariance Matrix
            R = cov(EB);
            Q = cov(EA');
            figure
            subplot(1,2,1)
            imagesc(Q);
            subplot(1,2,2)
            imagesc(R);
            %Final Plot
            figure
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
            figure
            for d = 1:D.Nd
                subplot(2,4,d)
                imagesc(D.Train{1,d}.FiringRate)
            end
        end
    end
end