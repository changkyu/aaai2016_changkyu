%  COV_D2F - 
function Sv = cov_d2f( Sv )

if isnumeric(Sv)
    Sv_in = Sv;
    
    Sv = cell(1,size(Sv,2));
    for i = 1:size(Sv,2)
        Sv{i} = diag(Sv_in(:,i));
    end
end

