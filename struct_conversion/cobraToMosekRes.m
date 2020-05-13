function res = cobraToMosekRes(res,keep,milp)

% build solution
if(res.time>=1e9)
    res.rcode='MSK_RES_TRM_MAX_TIME';
else
    res.rcode='crap';
end

if(milp)
    res.sol.int.solsta=res.stat;
    res.sol.int.prosta=res.stat;
    if(res.stat~=0)
        res.sol.int.xx=res.full(1:keep);
        %res.sol.int.xx=res.int;
        res.sol.int.pobjval=res.obj;
        res.sol.pobjval=res.obj;
    end
else
    res.sol.bas.solsta='FEASIBLE';
    res.sol.bas.prosta='OPTIMAL';
    if res.stat==0
        res.sol.bas.solsta='INFEASIBLE';
        res.sol.bas.prosta='INFEASIBLE';
    elseif res.stat==2
        res.sol.bas.prosta='UNBOUND';
    end
    if(res.stat~=0)
        res.sol.bas.xx=res.full(1:keep);
        res.sol.bas.pobjval=res.obj;
        res.sol.pobjval=res.obj;
        res.x=res.full(1:keep);
    end
end
end
