function [U, Sigma, V] = Bidiag_Francis_Step_Update_U_V( U, B, V )
    % size check for integrity. need B square and at least 2x2 for T
    % no instructions were given to make this work for non-square input
    [m, n] = size(B);
    if ((m ~= n) || ( m < 2 ))
        disp("Bidiag_Francis_Step_Update_U_V ERROR: input must be " + ...
            " square matrix of dim >= 2x2");
        Sigma = B;
        return;
    end

    % First Givens rotation is applied to T, but we only need T00, Tendend, T10
    x = [ B(1,1)^2 - (B(m-1,m)^2 + B(m,m)^2); ... % T00 - Tendend
          B(1,2) * B(1,1)]; % T10
    
    % Setup iteration for applying rotations
    doneBelow = false; % tracks bulge chasing below main diagonal
    doneAbove = false; % tracks bulge chasing above main diagonal
    bulgeBelowIdx = [1, 0]; % idx of bulge introduced below main diagonal. first will be 2,1, so init to 1,0
    bulgeAboveIdx = [0, 2]; % idx of bulge introduced above main diagonal. first will be 1,3, so init to 0,2
    xFirstIdx = [0, 0]; % will track idx of other element passed to Givens_rotation
    bulgeIdx = bulgeBelowIdx; % tracks where the bulge is for next iteration
    GIdx = [1, 2];     % to track where to put G in multiplying matrix
    applyRight = true; % we start by applying to the right
    iterCount = 0; % to track when we can update G indices
    Sigma = B;
    Vt = V';
    while ((doneBelow == false) || (doneAbove == false))
        % get this iteration's rotation
        if(iterCount > 0)
            x = [Sigma(xFirstIdx(1), xFirstIdx(2));...
                 Sigma(bulgeIdx(1),bulgeIdx(2))];
        end
        G = Givens_rotation(x);

        % Reference section 11.2.4
        GMatrix = eye(m);
        GMatrix(GIdx(1):GIdx(2), GIdx(1):GIdx(2)) = G;
        GMatrixT = GMatrix';

        % apply Gmatrix to chase bulge
        if (applyRight == true)
            % applying to the right introduces a bulge below main diagonal
            % A = U * B * V' => A = U * B * G * G' * V'
            % Update Sigma: optimize by only multiplying necessary entries
            % protect index overflow with if...else
            if ((bulgeIdx(1) + 2) <= m)
                Sigma(bulgeIdx(1):bulgeIdx(1)+2, GIdx(1):GIdx(2)) = ...
                    Sigma(bulgeIdx(1):bulgeIdx(1)+2, GIdx(1):GIdx(2)) * G;
            else
                GMatrix = eye(m);
                GMatrix(GIdx(1):GIdx(2), GIdx(1):GIdx(2)) = G;
                Sigma = Sigma * GMatrix;
            end
            
            % Update V': optimize by only getting necessary row updated
            Vt(GIdx(1):GIdx(2), :) = GMatrixT(GIdx(1):GIdx(2), :) * Vt;
            bulgeBelowIdx = bulgeBelowIdx + [1, 1];

            % apply left next time
            bulgeIdx = bulgeBelowIdx;
            xFirstIdx = bulgeIdx - [1, 0];
            applyRight = false; 

            % check if we are done chasing below the matrix 
            if(bulgeBelowIdx(1) >= m)
                doneBelow = true;
            end
        else
            % applying to the left introduces a bulge above main diagonal
            % A = U * B * V' => A = U * G * G'* B * V'
            % Update Sigma: optimize by only multiplying necessary entries
            if((bulgeIdx(2) <= (m-2)))
                Sigma(GIdx(1):GIdx(2), bulgeIdx(2):bulgeIdx(2)+2) = ...
                    G' * Sigma(GIdx(1):GIdx(2), bulgeIdx(2):bulgeIdx(2)+2);
            else % last case only touches the bottom right 2x2 of Bi_next
                Sigma(GIdx(1):GIdx(2), (m-1):m) = ...
                    G' * Sigma(GIdx(1):GIdx(2), (m-1):m);
            end

            % Update U: optimize by only getting necessary columns updated
            U(:, GIdx(1):GIdx(2)) = U * GMatrix(:, GIdx(1):GIdx(2));
            bulgeAboveIdx = bulgeAboveIdx + [1, 1];

            % apply right next time
            bulgeIdx = bulgeAboveIdx;
            xFirstIdx = bulgeIdx - [0, 1];
            applyRight = true;

            % check if we are done chasing above the matrix 
            if(bulgeAboveIdx(2) > m)
                doneAbove = true;
            end
        end

        % update G indices and x for next iteration if we aren't done
        iterCount = iterCount + 1;
        if(mod(iterCount,2) == 0)
                GIdx = GIdx + [1,1];
        end
    end

    % Get to correct output format
    V = Vt';
end