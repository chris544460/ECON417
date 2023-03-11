using LinearAlgebra
include("io2.jl")
include("creategrid1.jl")
include("spline2.jl")

function interp(x,yspline;calcy=true,calcyp=false,calcydp=false)

    one = 1.0
    
    npts = yspline.n
 
    klo = 1
    khi = npts
    while ((khi-klo) > 1)
       k = Int(trunc((khi+klo)/2))
       if (yspline.x[k] > x)
          khi = k
       else
          klo = k
       end
    end
    
    h = yspline.x[khi] - yspline.x[klo]
    a = (yspline.x[khi] - x)/h
    b = (x - yspline.x[klo])/h
    asq = a*a
    bsq = b*b
 
    if calcy
       y = a*yspline.y[klo] + b*yspline.y[khi] + 
           ((asq*a-a)*yspline.ydp[klo] 
            +(bsq*b-b)*yspline.ydp[khi])*(h*h)/6
    else
       y = 0.0
    end           
 
    if calcyp
       yp = (yspline.y[khi]-yspline.y[klo])/h - 
            (3*asq-one)/6*h*yspline.ydp[klo] + 
            (3*bsq-one)/6*h*yspline.ydp[khi]
    else
       yp = 0.0
    end           
 
    if calcydp
       ydp = a*yspline.ydp[klo] + b*yspline.ydp[khi]
    else
       ydp = 0.0
    end
    
    interp = (y,yp,ydp)      
       
 end 

# define interp1d function to be used later 
function interp1d(x,yspline;calcy=true,calcyp=false,calcydp=false)
    return interp(x,yspline;calcy=calcy,calcyp=calcyp,calcydp=calcydp)
end 
 


 # brent is a function that finds the minimum of a function f(x) on the interval
function brent(ax,bx,cx,f,tol,zeps)

    itmax = 100
    cgold = 0.3819660
    toler = 1.0e-6
    zero = 0.0
 
    a = min(ax,cx)
    b = max(ax,cx)
    v = bx
    w = v
    x = v
    e = zero
    fx = f(x)
    fv = fx
    fw = fx
    
    d = 0.0
    
    for iter in 1:itmax
       xm = 0.5*(a+b)
       tol1 = tol*abs(x) + zeps
       tol2 = 2.0*tol1
       if (abs(x-xm) <= (tol2-0.5*(b-a)))
          fmin = fx
          xmin = x
          niter = iter
          return fmin,xmin,niter
       end   
       if (abs(e) > tol1) 
          r = (x-w)*(fx-fv)
          q = (x-v)*(fx-fw)
          p = (x-v)*q - (x-w)*r
          q = 2.0*(q-r)
          if (q > zero) 
             p = -p
          end   
          q = abs(q)
          etemp = e
          e = d
          if ( (abs(p) >= abs(0.5*q*etemp)) |
               (p <= (q*(a-x))) | (p >= (q*(b-x))) ) 
             if (x >= xm) 
                e = a - x
             else
                e = b - x
             end
             d = cgold*e
          else   
             d = p/q
             u = x + d
             if ( ((u-a) < tol2) | ((b-u) < tol2) )
                d = abs(tol1)*sign(xm-x)
             end   
          end  
       else
          if (x >= xm) 
             e = a - x
          else
             e = b - x
          end      
          d = cgold*e
       end
       if (abs(d) >= tol1) 
          u = x + d
       else
          u = x + abs(tol1)*sign(d)
       end
       fu = f(u)
       if (fu <= fx) 
          if (u >= x) 
             a = x
          else
             b = x
          end
          v = w
          fv = fw
          w = x
          fw = fx
          x = u
          fx = fu
       else
          if (u < x) 
             a = u
          else
             b = u
          end
          if ( (fu <= fw) | (abs(w-x) <= toler) )
             v = w
             fv = fw
             w = u
             fw = fu
          elseif ( (fu <= fv) | (abs(v-x) <= toler) | 
                   (abs(v-w) <= toler) ) 
             v = u
             fv = fu
          end
       end 
    end
 
    wait("Maximum number of iterations exceeded in brent.")
 
    fmin = fx
    xmin = x
    niter = iter
 
    return fmin,xmin,niter
 
 end
 

# Define the production function
function f(k)
    return k^alpha # Cobb Douglas production function
end

# Define the utility optimal_value_function
function U(c)
    return log(c) # log utility function
end

# Define the optiml decision rule
function q(k, yspline, k_min, k_max)
    q = brent(k_min, k_max, k_max, x -> U(f(x) - x) + yspline(x), 1e-6, 1e-6)[2]
    return q
end

# Define the function to compute the optimal value function and the optimal decision rule

function optimal_value_function(alpha, beta, delta, N, tolerance)
    k_ss = (alpha / (1 / beta - 1 + delta))^(1 / (1 - alpha)) # steady-state capital stock
    k = range(0.5 * k_ss, 1.5 * k_ss, length = N) # grid for capital
    v = zeros(N) # initial set of values for the value function (step i)
    while true
        v_new = zeros(N) # new set of values for the value function
        yspline = makespline(k, v) # create a cubic spline of the value function
    
        for i in 1:N
            ki = k[i] # ki is the level of capital at the i-th grid point
            k_min = ki^alpha / (1 + delta)
            k_max = f(ki) + (1 - delta) * ki
            k_star = q(ki, yspline, k_min, k_max)
            v_new[i] = -U(f(ki) - k_star + (1 - delta) * ki) - beta * interp1d(k_min, k_max, yspline, f(ki)) # compute the new value function (step ii)
        end
    
        distance = maximum(abs.(v_new - v)) # compute the distance between the updated value function and the previous value function (step iv)
        v = v_new # update the value function (step iv)
    
        if distance < tolerance
            break
        end
    end
    
    g_return = [q(k[i], yspline, k[i]^alpha / (1 + delta), f(k[i]) + (1 - delta) * k[i]) for i in 1:N] # optimal choice

    return v, g_return
end

# Set parameter values
alpha = 0.4
beta = 0.9
delta = 0.1
tolerance = 1e-6

# Compute the steady-state capital stock
k_star = ((1 / beta - 1 + delta) / alpha)^(1 / (alpha - 1)) # steady-state capital stock

# Set up the grid for capital
N = 31
k_min = 0.5 * k_star # minimum value of capital
k_max = 1.5 * k_star # maximum value of capital
k = range(k_min, k_max, length = N) # grid for capital

# Compute the optimal value function and the optimal decision rule
v, g = optimal_value_function(alpha, beta, delta, N, tolerance)

# Plot the optimal value function
plot(k, v, color = "blue", linewidth = 2, label = "optimal value function")
# plot(k, U(f(k) - g + (1 - delta) * k), color = "red", linewidth = 2, label = "optimal value function")
    
    