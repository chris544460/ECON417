# Code for generate Gauss-Hermite weights and abscissas.

include("io2.jl")

# This function computes the weights, w, and abscissas, x, for n-point Gauss-Hermite quadrature.
# The function returns the two arrays w and x, both of size n.
function gauher(n)

   maxit = 10
   eps = 3.0e-14
   pim4 = 0.7511255444649425
   pi = 3.14159265358979
   
   pp = 0.0
   z = 0.0
   
   x = zeros(n)
   w = zeros(n)

   m = trunc(Int,(n+1)/2)
   
   for i in 1:m

      if (i == 1) 
         xn = 2*n + 1
         z = sqrt(xn) - 1.85575*(xn^(-0.16667))
      elseif (i == 2) 
         z = z - 1.14*(n^0.426)/z
      elseif (i == 3) 
         z = 1.86*z - 0.86*x[1]
      elseif (i == 4)
         z = 1.91*z - 0.91*x[2]
      else
         z = 2.0*z - x[i-2]
      end
      
      done = false
      
      for its in 1:maxit
         iter = its
         p1 = pim4
         p2 = 0.0
         for j in 1:n
            p3 = p2
            p2 = p1
            p1 = z*sqrt(2.0/j)*p2 - sqrt((j-1)/j)*p3
         end
         pp = sqrt(2.0*n)*p2
         z1 = z
         z = z1 - p1/pp
         if (abs(z-z1) <= eps) 
            done = true
            break
         end   
      end
      
      if !done
         wait("Too many iterations in gauher")
      else   
         x[i] = z
         x[n+1-i] = -z
         w[i] = 2.0/(pp*pp)
         w[n+1-i] = w[i]
      end   
      
   end   

   return (w,x)

end 
