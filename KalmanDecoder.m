classdef KalmanDecoder
    
    properties 
    %Model
    Nx = 1;
    A = [];
    Q = [];
    H = [];
    R = [];
    P = [];
    C =  [];
    C1 = [];
    C2 = [];
    C3 = [];
    
    end
    
    methods
        %% Constructor
        function Filter = KalmanDecoder(A,Q,H,R)
            %Model
            Filter.Nx = size(A);
            Filter.A = A; 
            Filter.Q = Q; 
            %Observation
            Filter.H = H; 
            Filter.R = R;
            %Corretion
            [n,~] = size(Q);
            Filter.P = zeros(n,n);
        end
        %% Update
        function [K,XNew] = KalmanUpdate(K,XOld,Z)
            %Estimation
            XEst = K.A*XOld;
            %Update System Variance
            Pm = K.A*K.P*K.A' + K.Q;
            %Kalman Gain
            Kg = Pm*K.H'*(K.H*Pm*K.H' + K.R)^-1;
            %Correction
            XNew = XEst + Kg*(Z - K.H*XEst);
            K.C =  [K.C,XEst];
            K.C1 = [K.C1,Kg*(Z - K.H*XEst)];
            K.C2 = [K.C2,Kg*(Z)];
            K.C3 = [K.C3,Kg*(K.H*XEst)];
            
            %Update Error Covariance
            K.P  = (eye(K.Nx)-Kg*K.H)*Pm;
            
        end
        function [Filter] = ReInitialise(Filter)
            [n,~] = size(Filter.P);
            Filter.P = eye(n);
        end
    end
    
    methods (Static)
        function [H] = Observation(B)
            %Transition
            T = [1 0 0 0 0 0
                 0 0 0 1 0 0
                 0 0 0 0 1 0
                 0 0 0 0 0 1];
            %Observation
            C = [1         0         0         0;
                 0.4643    0.3827    0.1447   -0.0379;
                -0.5859    0.0887    0.4535    0.0478;
                 1.6311   -0.0581    0.0483    0.4165;];
                 
            %C = eye(4);
            %B = B*C;
            H = B*T;
        end
        function [H,Zo] = PolarObservation(B,r)
            %Prepare
            T = [1 0 -1  0;
                 0 r  0 -r];
            H = B(:,2:3)*T;
            Zo = B(:,1);
        end
    end
end