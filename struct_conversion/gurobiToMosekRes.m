function res = gurobiToMosekRes(res,keep,milp)
res.rcode=1000;

if(res.runtime>=1e9)
    res.rcodestr='MSK_RES_TRM_MAX_TIME';
else
    res.rcodestr='crap';
end
try
    if(milp)
        res.sol.int.solsta=res.status;
        res.sol.int.prosta=res.status;
        if(~strcmp(res.status,'INFEASIBLE'))
            res.sol.int.xx=res.x(1:keep);
            res.sol.int.pobjval=res.objval;
            res.sol.pobjval=res.objval;
        end
    else
        res.sol.bas.solsta=res.status;
        res.sol.bas.prosta=res.status;
        if(~strcmp(res.status,'INFEASIBLE'))
            res.sol.bas.xx=res.x(1:keep);
            res.sol.bas.pobjval=res.objval;
            res.sol.pobjval=res.objval;
        end
    end
    
    res.x=res.x(1:keep);
catch
    res.rcode=0;
    res.rmsg='';
    res.rcodestr='MSK_RES_OK';
    res.sol.bas.prosta='PRIMAL_INFEASIBLE';
end
end
