function findLastWorking(testClass,funName)
    % Use to find last working revision, e.g.,
    % findLastWorking(HypertoroidalFourierDistributionTest,'testMultiply2D');
    [~,output]=system('git rev-list --branches master --pretty=oneline | cut -c1-40');
    allVerisions=regexp(output,'\n','split');
    for currVers=allVerisions
        system(['git checkout ',currVers{:}]);
        result=run(testClass,funName);
        if ~result.Failed
            disp(['Working with version ',currVers{:}]);
            return
        end
    end
end