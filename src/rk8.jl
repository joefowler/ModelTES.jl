# For the system of ODEs du/dx = df(x,u), take one step with eighth order
# (ninth order locally) Runge Kutta.  y contains initial values at the
# start point x=t.  The final values at x=t+h are given in ys.  work
# is scratch space of size 7n, where n is the dimension of y and ys.
# Method of G. J. Cooper and J. H. Verner, SIAM Journal of Numerical
# Analysis 9, pp. 389-405 (1972).
function rk8!(df, t, h, y, ys, work)
    n = length(y)
    rk8aux(df, t, h, y, ys, view(work,1:n),
           view(work,  n+1:2*n), view(work,2*n+1:3*n), view(work,3*n+1:4*n),
           view(work,4*n+1:5*n), view(work,5*n+1:6*n), view(work,6*n+1:7*n))
end
function rk8aux(df, t, h, y, ys, yp, f1,f2,f3,f4,f5,f6)
    n = length(y)
    df(t, y, yp)
    for k = 1:n
        ys[k] = y[k] + h * (1/2) * yp[k]
    end
    df(t+h*(1/2), ys, f1)
    for k = 1:n
        ys[k] = y[k] + h * (1/4) * (yp[k] + f1[k])
    end
    df(t+h*(1/2), ys, f2)
    for k = 1:n
        ys[k] = y[k] + h * ((1/7) * yp[k] +
                            (-7+3*sqrt(21))/98 * f1[k] +
                            (21-5*sqrt(21))/49 * f2[k] )
    end
    df(t+h*(7-sqrt(21))/14, ys, f1)
    for k = 1:n
        ys[k] = y[k] + h * ((11-sqrt(21))/84 * yp[k] +
                            (18-4*sqrt(21))/63 * f2[k] +
                            (21+sqrt(21))/252 * f1[k])
    end
    df(t+h*(7-sqrt(21))/14, ys, f3)
    for k = 1:n
        ys[k] = y[k] + h * ((5-sqrt(21))/48 * yp[k] +
                            (9-sqrt(21))/36 * f2[k] +
                            (-231-14*sqrt(21))/360 * f1[k] +
                            (63+7*sqrt(21))/80 * f3[k])
    end
    df(t+h*(1/2), ys, f4)
    for k = 1:n
        ys[k] = y[k] + h * ((10+sqrt(21))/42 * yp[k] +
                            (-432-92*sqrt(21))/315 * f2[k] +
                            (633+145*sqrt(21))/90 * f1[k] +
                            (-504-115*sqrt(21))/70 * f3[k] +
                            (63+13*sqrt(21))/35 * f4[k])
    end
    df(t+h*(7+sqrt(21))/14, ys, f2)
    for k = 1:n
        ys[k] = y[k] + h * ((1/14) * yp[k] +
                            (14+3*sqrt(21))/126 * f3[k] +
                            (13+3*sqrt(21))/63 * f4[k] +
                            (1/9) * f2[k])
    end
    df(t+h*(7+sqrt(21))/14, ys, f1)
    for k = 1:n
        ys[k] = y[k] + h * ((1/32) * yp[k] +
                            (91+21*sqrt(21))/576 * f3[k] +
                            (11/72) * f4[k] +
                            (-385+75*sqrt(21))/1152 * f2[k] +
                            (63-13*sqrt(21))/128 * f1[k])
    end
    df(t+h*(1/2), ys, f5)
    for k = 1:n
        ys[k] = y[k] + h * ((1/14) * yp[k] +
                            (1/9) * f3[k] +
                            (-733+147*sqrt(21))/2205 * f4[k] +
                            (515-111*sqrt(21))/504 * f2[k] +
                            (-51+11*sqrt(21))/56 * f1[k] +
                            (132-28*sqrt(21))/245 * f5[k])
    end
    df(t+h*(7-sqrt(21))/14, ys, f6)
    for k = 1:n
        ys[k] = y[k] + h * ((-42-7*sqrt(21))/18 * f3[k] +
                            (-18-28*sqrt(21))/45 * f4[k] +
                            (-273+53*sqrt(21))/72 * f2[k] +
                            (301-53*sqrt(21))/72 * f1[k] +
                            (28+28*sqrt(21))/45 * f5[k] +
                            (49+7*sqrt(21))/18 * f6[k])
    end
    df(t+h, ys, f2)
    # approximate solution at t+h
    for k = 1:n
        ys[k] = y[k] + h * ((1/20) * yp[k] +
                            (49/180) * f1[k] +
                            (16/45) * f5[k] +
                            (49/180) * f6[k] +
                            (1/20) * f2[k])
    end
end
