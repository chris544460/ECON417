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



function optimal_value_function(N, k_steadyState, F, U; alpha=0.4, beta=0.9, delta=0.1, toler=1.0e-10, maxiter=5000)

    v0 = zeros(N)
    A = creategrid(0.5*k_steadyState, 1.5*k_steadyState, N)
    spline_ = makespline(A, v0)
    iter = 1

    while iter < maxiter
        iter += 1
        V1 = zeros(N)
        kstar = zeros(N)

        for i in 1:N
            fk = F(A[i])
            (v_min, k_min) = brent(0, fk / 2, fk, kprime -> U(fk - kprime + 1e-10) + beta * interp(kprime, spline_)[1], 1.0e-8, 1.0e-8)

            V1[i] = max(v_min, 0)
            kstar[i] = k_min

        end

        distance = maximum(abs.(V1 - v0))

        if distance < toler
            return V1, kstar, spline_
        else
            v0 = V1
            spline_ = makespline(A, v0)
        end

    end

    return "Value function iteration did not converge"

end


