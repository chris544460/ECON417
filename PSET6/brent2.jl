include("io2.jl")

# This code, drawn from Numerical Recipes, implements Brent's method (without derivative) to
# minimize the function f.  The arguments ax, bx, and cx must satisfy: ax < bx < cx, f(ax) > f(bx),
# and f(cx) > f(bx).  The function returns (fmin,xmin,niter), where fmin is the minimized value of
# the function f, xmin is the minimand, and niter is the number of iterations required to
# find the minimum.  When using brent, good choices for the arguments tol and zeps are:
# tol = zeps = 1.0e-8.  Here is an example of how to "call" brent:
#
# fmin,xmin,niter = brent(ax,bx,cx,f,1.0e-8,1.0e-8).


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
