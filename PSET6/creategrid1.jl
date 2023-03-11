# Create an equally-spaced grid with npts points on the interval [xlow,xhigh].

function creategrid(xlow,xhigh,npts)

   xinc = (xhigh - xlow)/(npts-1)

   return collect(xlow:xinc:xhigh)

end
