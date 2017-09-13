
module PsychrometricCharts

# package code goes here

using Psychro

function calcpoints(Pa=101325.0, Tlim=[0.0, 50.0], nT=100,
                    Wlim=[0.0, 30.0], nW=100)

    Tmax = maximum(Tlim) + 273.15
    Tmin = minimum(Tlim) + 273.15
    dT = (Tmax - Tmin) / (nT-1)
    T = Tmin:dT:Tmax

    
    Wmax = maximum(Wlim)/1000
    Wmin = minimum(Wlim)/1000
    dW = (Wmax - Wmin) / (nW-1)
    


    h = zeros(nW, nT)
    s = zeros(nW, nT)
    r = zeros(nW, nT)
    Td = zeros(nW, nT)
    w = zeros(nW, nT)
    q = zeros(nW, nT)
    ρ = zeros(nW, nT)
    Tb = zeros(nW, nT)
    xv = zeros(nW, nT)
    q = zeros(nW, nT)
    Z = zeros(nW, nT)


    wsat = humrat.(MoistAir, T, RelHum, 1.0, Pa)
    
    for i = 1:nT
        w1 = linspace(Wmin, max(wsat[i], Wmax), nW)
        for k = 1:nW
            h[i,k] = enthalpy(MoistAir, T[i], HumRat, w1[k], Pa)
            s[i,k] = entropy(MoistAir, T[i], HumRat, w1[k], Pa)
            r[i,k] = relhum(MoistAir, T[i], HumRat, w1[k], Pa)
            
        end
    end
end

function pretty(x, n=8, dx=[1.0, 2.0, 5.0])
    xmin = minimum(x)
    xmax = maximum(x)
    xr = xmax-xmin
    dx0 = xr / n
    logdx = [log10(y) - floor(log10(y))  for y in dx]
    pwr = floor(log10(dx0))
    l0 = log10(dx0) - pwr
    
    val, idx = findmin(abs.(abs(l0) .- logdx))
    #return val, idx, dx[idx], (l0, pwr)
    dx1 = dx[idx]*10^pwr
    
    xlim1 = round(xmin/dx1, RoundDown)*dx1
    xlim2 = round(xmax/dx1, RoundUp)*dx1
    
    return xlim1:dx1:xlim2
    
end
   
function drawlines(humfun, HumType, n=10, Tc, W, Pa=101325.0)

    Tmin = minimum(Tc) + 273.15
    Tmax = maximum(Tc) + 273.15
    Wmin = minimum(W)
    Wmax = maximum(W)

    wsat1 = humrat(MoistAir, Tmin, RelHum, 1.0, Pa)
    wsat2 = humrat(MoistAir, Tmax, RelHum, 1.0, Pa)

    wmin = min(Wmin, wsat1)
    wmax = min(wsat2, Wmax)

    p1 = humfun(MoistAir, Tmin, HumRat, wmin, Pa)
    p2 = humfun(MoistAir, Tmin, HumRat, wmax, Pa)
    p3 = humfun(MoistAir, Tmax, HumRat, wmin, Pa)
    p4 = humfun(MoistAir, Tmax, HumRat, wmax, Pa)

    pmin = min(p1, p2, p3, p4)
    pmax = max(p1, p2, p3, p4)

    
    
    

end
function relhumlines(Tc, r, Pa)
    nT = size(T,1)
    nr = size(r,1)

    w = zeros(nT, nr)
    h = zeros(nT, nr)

    T0 = 273.15

    for i in 1:nr
        for k in 1:nT
            w[k,i] = humrat(MoistAir, Tc[k] + T0, RelHum, r[i], Pa)
            h[k,i] = enthalpy(MoistAir, T[k] + T0, RelHum, r[i], Pa)
        end
    end

    return w, h
end

    


end # module

            

