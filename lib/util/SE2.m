%% SE2 Group of proper rigid transformations in the plane.
%
% This class implements the group of proper rigid transformations and
% its different representations. Internally, a proper rigid
% transformation is stored using a rotation and a subsequent
% translation.

classdef SE2   
    properties
        alpha double % Angle of rotation.
        translation double % Translation (which is performed after the rotation)
    end
    
    methods        
        function this = SE2(alpha, translation)
            % Constructor of the class.
            %
            %   Parameters:
            %       alpha - rotation (as a radian)
            arguments
                alpha (1,1) double
                translation (2,1) double
            end
            this.alpha = alpha;
            this.translation = reshape(translation,2,1);
        end
        
        function dq = asDualQuaternion(this, style)
            % Returns the dual quaternion representation.
            %
            %   Returns the dual quaternion representation of the current
            %   proper rigid transformation.
            %
            %   Parameters:
            %       style - 'full' represents a dual quaternion as an 8 dim
            %               vector. 'compact' (default) represents it as an 4
            %               dim vector.
            arguments
                this (1,1) SE2
                style char = 'compact'
            end
            
            switch style
                case 'full'
                    dq = zeros(8,1);
                    dq(1) = cos(this.alpha/2);
                    dq(4) = sin(this.alpha/2);
                    dq(6) = 0.5*(cos(this.alpha/2) * this.translation(1) ...
                        + sin(this.alpha/2)*this.translation(2));
                    dq(7) = 0.5*(-sin(this.alpha/2) * this.translation(1) ...
                        + cos(this.alpha/2)*this.translation(2));
                otherwise % compact style is default.
                    dq = zeros(4,1);
                    dq(1) = cos(this.alpha/2);
                    dq(2) = sin(this.alpha/2);
                    dq(3) = 0.5*(cos(this.alpha/2) * this.translation(1) ...
                        + sin(this.alpha/2)*this.translation(2));
                    dq(4) = 0.5*(-sin(this.alpha/2) * this.translation(1) ...
                        + cos(this.alpha/2)*this.translation(2));
            end
        end
        
        function dqMat = asDualQuaternionMatrix(this)
            % Returns the dual quaternion Matrix representation.
            %
            % Returns the dual quaternion representation of the current
            % proper rigid transformation. A compact 4 x 4 representation is
            % used.
            
            dq = this.asDualQuaternion();
            dqMat = SE2.dualQuaternionToMatrix(dq);
        end
        
        function hMat = asHomogenousMatrix(this)
            % Represents transformation as Homogenous matrix.
            hMat = [this.rotationPartAsRotationMatrix this.translation;
                0 0 1];            
        end
        
        function rotMat = rotationPartAsRotationMatrix(this)
            % Returns a rotation matrix.
            rotMat = [cos(this.alpha) -sin(this.alpha);
                sin(this.alpha)  cos(this.alpha)];
        end       
    end
    
    methods (Static)
        function se2 = fromDualQuaternion(dq)
            arguments
                dq (4,1) double
            end
            [angle,pos] = AbstractSE2Distribution.dualQuaternionToAnglePos(dq);
            se2 = SE2(angle,pos);
        end
        function res = dualQuaternionMultiply(dq1, dq2)
            % Multiplies two dual quaternions.
            assert(isequal(size(dq1),[4 1]));
            assert(isequal(size(dq1),[4 1]));
            
            res = SE2.matrixToDualQuaternion(...
                SE2.dualQuaternionToMatrix( dq1 ) ...
                * SE2.dualQuaternionToMatrix( dq2 ) );
        end
        
        function dqMat = dualQuaternionToMatrix(dq)
            % Transforms a dual quaternion into a matrix.
            %
            %   A compact Matrix representation of the dual quaternion is
            %   returned.          
            assert(all(size(dq) == [4 1]));
            
            dqMat = [ dq(1)  dq(2)    0      0
                -dq(2)  dq(1)    0      0
                -dq(3)  dq(4) dq(1) -dq(2)
                -dq(4) -dq(3) dq(2)  dq(1) ];
        end
        
        function dq = matrixToDualQuaternion(dqMat)
            % Transform our matrix representation to dq.
            %
            % A compact Matrix representation of the dual quaternion is
            % returned.
            assert(isequal(size(dqMat),[4 4]));
            
            dq=zeros(4,1);
            
            dq(1) =  dqMat(1,1);
            dq(2) = -dqMat(2,1);
            dq(3) = -dqMat(3,1);
            dq(4) = -dqMat(4,1);
        end 
    end
end

