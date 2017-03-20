classdef KalmanDecoder
    
    properties 
    %Model
    Nx = 1;
    A = [];
    Q = [];
    H = [];
    R = [];
    P = [];
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
            Filter.P = eye(n);
        end
        %% Update
        function [Filter,XNew] = KalmanUpdate(Filter,XOld,Z)
            %Estimation
            XEst = Filter.A*XOld;
            %Update System Variance
            Pm = Filter.A*Filter.P*Filter.A' + Filter.Q;
            %Kalman Gain
            Kg = Pm*Filter.H'*(Filter.H*Pm*Filter.H' + Filter.R)^-1;
            %Correction
            XNew = XEst + Kg*(Z - Filter.H*XEst);
            %Update Error Covariance
            Filter.P  = (eye(Filter.Nx)-Kg*Filter.H)*Pm;
        end
        function [Filter] = ReInitialise(Filter)
            [n,~] = size(Filter.P);
            Filter.P = eye(n);
        end
    end
    
    methods (Static)
        function [H,Zo] = Observation(B)
            %Prepare
            T = [1 0 0 -1  0  0;
                 0 1 0  0 -1  0
                 0 0 1  0  0  0];
            H = B(:,2:4)*T;
            Zo = B(:,1);
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