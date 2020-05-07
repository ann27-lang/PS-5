# CHEME 5440/7770
# Problem 2
include("PhasePortraitV2.jl")

function odefun(ca, cr)
    # Parameters
    ror = 1.0
    rr = 100.0
    roa = 100.0
    ra = 5000.0
    da = 30.0

    u = @. -da*ca + (roa+ra*ca^2)/(1+ca^2+cr^2)
    v = @. -cr + (ror+rr*ca^2)/(1+ca^2)

    return(u,v)
end


tspan = (0.0, 50.0)

xlimca = (0, 200, 50)
ylimcr = (0, 100, 50)
xo = ([1.0, 10.0],)

phaseplot(odefun,xlimca,ylimcr,clines=true,xinit=xo,t=tspan,norm=true,scale=0.5)
