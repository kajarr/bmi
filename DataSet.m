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
        function D = DataSet(data,Ntr)  
            %Look for every direction
            for d = 1:D.Nd
                %Temp Struct Time
                temp = struct;
                S = zeros(D.Nn,D.Nt-1);
                X = zeros(2,D.Nt);
                V = zeros(2,D.Nt-1);  
                %Look for every trial
                for n = Ntr(1):Ntr(2)
                    %Local Data    
                    trial = data(n,d);
                    s = trial.spikes(:,2:D.Nt);
                    x = trial.handPos(1:2,1:D.Nt);
                    v = diff(x')';
                    %Spike, Position, Velocity,TrialMean
                    %Trial
                    S = S + s;
                    X = X + x;
                    V = V + v;
                end
            %Average   
            N = Ntr(2) - Ntr(1) + 1;   
            temp.Spikes = S/N;
            temp.Position = X/N;
            temp.Velocity = V/N;
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
        function D = UnitVector(D,Field)
            for d = 1:D.Nd
                for t = 1:D.Nt-1
                V = D.Dir{d}.(Field)(:,t);
                v = norm(V);
                V = V/v;
                D.Dir{d}.(Field)(:,t) = V;
                end
            end
        end
        %% Output
        function [W,E] = GetPreDirection(D,a)
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
            B = (V'*V+a(1)*ones(2,2))^-1*V'*F;
            E = F-V*B;
            %Variance of the Residuals
            s = ones(1,D.Nn)./var(E);
            S = diag(s);
            W = pinv(B*S*B' + a(2)*ones(2,2))*B*S;
        end
        
        %% Test
        function [] = MeanTest(W,Test)
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
                X = zeros(2,Test.Nt);
                V = zeros(2,Test.Nt-1);
                    for t = 1:Test.Nt-1
                        %Regress
                        V(:,t) = W*(F(:,t));
                        X(:,t+1) = X(:,t) + V(:,t);    
                    end
                figure(h1) 
                plot(X(1,:),X(2,:))
                
                figure(h2)
                subplot(2,4,d)
                hold on
                l = plot(V');
                set(l(1),'Color','b')
                set(l(2),'Color','b','LineStyle','--')
                l = plot(Test.Dir{d}.Velocity');
                set(l(1),'Color','r')
                set(l(2),'Color','r','LineStyle','--')
            end 
        end
    end
end