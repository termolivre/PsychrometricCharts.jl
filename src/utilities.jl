


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


function solvezero(f, x0, x1, xorder=0.0, EPS=1e-8, MAXITER=100)

    if xorder==0.0
        xorder = 1.0
    end
    y0 = f(x0)
    y1 = f(x1)
    x = 0.0
    for i = 1:MAXITER
        a = (y1-y0) / (x1-x0)
        b = y0 - a*x0
        x = -b/a
        if abs(x-x1) < EPS*abs(xorder)
            return x
        end
        x0 = x1
        y0 = y1
        x1 = x
        y1 = f(x1)
    end

    error("solvezero did not converge")
end
