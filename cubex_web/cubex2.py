#!/usr/bin/python
import cgi,os,time
import math

form = cgi.FieldStorage()

def cubex(inputarray):
    result = {}
    n1111 = float(inputarray[0])
    n1112 = float(inputarray[1])
    n1122 = float(inputarray[2])
    n1211 = float(inputarray[3])
    n1212 = float(inputarray[4])
    n1222 = float(inputarray[5])
    n2211 = float(inputarray[6])
    n2212 = float(inputarray[7])
    n2222 = float(inputarray[8])
    n = (n1111 + n1112 + n1122 + n1211 + n1212 + n1222 + n2211 + n2212 + n2222)
    result["n"] = n
    p = ((n1111+n1112+n1122)*2.0+(n1211+n1212+n1222))/(2.0 * n)
    q = ((n1111+n1211+n2211)*2.0+(n1112+n1212+n2212))/(2.0 * n)
    n11 = (2.0*n1111 + n1112 + n1211) 
    n12 = (2.0*n1122 + n1112 + n1222)
    n21 = (2.0*n2211 + n2212 + n1211)
    n22 = (2.0*n2222 + n2212 + n1222) 
    a0 = -n11*p*q
    a1 = -n11*(1.0 - 2.0*p - 2.0*q) - n1212*(1.0 - p - q) + 2.0*n*p*q
    a2 = 2.0*n*(1.0 - 2.0*p - 2.0*q) - 2.0*n11 - n1212
    a3 = 4.0 * n
    minhap = n11 / (2.0 * float(n))
    maxhap = (n11 + n1212) / (2.0 * float(n))
    result["minhap"] = minhap
    result["maxhap"] = maxhap
    if p < 1.0 and p > 0.0:
        result["hwchisnp1"] = ((((n1111 + n1112 + n1122) - ((p ** 2)*n)))**2)/((p ** 2)*n)+\
                    ((((n1211 + n1212 + n1222) - ((2 * p * (1.0-p))*n)))**2)/((2 * p * (1.0-p))*n)+\
                    ((((n2211 + n2212 + n2222) - (((1-p) ** 2)*n)))**2)/(((1-p) ** 2)*n)
    else:
        result["hwchisnp1"] = 0.0
    if q < 1.0 and q > 0.0:
        result["hwchisnp2"] = ((((n1111 + n1211 + n2211) - ((q ** 2)*n)))**2)/((q ** 2)*n)+\
                    ((((n1112 + n1212 + n2212) - ((2 * q * (1.0-q))*n)))**2)/((2 * q * (1.0-q))*n)+\
                    ((((n1122 + n1222 + n2222) - (((1-q) ** 2)*n)))**2)/(((1-q) ** 2)*n)
    else:
        result["hwchisnp2"] = 0.0
    
    a = a3
    b = a2
    c = a1
    dee = a0
    
    xN = -b/(3.0*a)
    d2 = (math.pow(b,2)-3.0*a*c)/(9*math.pow(a,2))
    yN = a * math.pow(xN,3) + b * math.pow(xN,2) + c * xN + dee
    yN2 = math.pow(yN,2)

    h2 = 4 * math.pow(a,2) * math.pow(d2,3)
    result["realnoproblem"] = 0
    if abs(yN2-h2) <= 0.0000001:
        result["realnoproblem"] = 1

    if yN2 > h2 and result["realnoproblem"] == 0:
        # option 1
        number1 = 0.0
        number2 = 0.0
        if (1.0/(2.0*a)*(-yN + math.pow((yN2 - h2),0.5))) < 0:
            number1 = -math.pow(-(1.0/(2.0*a)*(-yN + math.pow((yN2 - h2),0.5))),1.0/3.0)
        else: number1 = math.pow((1.0/(2.0*a)*(-yN + math.pow((yN2 - h2),0.5))),1.0/3.0)
        if (1.0/(2.0*a)*(-yN - math.pow((yN2 - h2),0.5))) < 0:
            number2 = -math.pow(-(1.0/(2.0*a)*(-yN - math.pow((yN2 - h2),0.5))),1.0/3.0)
        else: number2 = math.pow((1.0/(2.0*a)*(-yN - math.pow((yN2 - h2),0.5))),1.0/3.0)
        alpha = xN + number1 + number2
        result["inputdata"] = str(inputarray)
        result["alpha"] = alpha
        result["beta"] = "Not a real root"
        result["gamma"] = "Not a real root"
        result["p"] = p
        result["q"] = q
        result["pq"] = str(p * q)
        if result["alpha"] >= minhap - 0.00001 and result["alpha"] <= maxhap + 0.00001:
            result["alphaposs"] = 1
            f11 = result["alpha"]
            f12 = p - f11
            f21 = q - f11
            f22 = 1 - (f11 + f12 + f21)
            D = (f11 * f22) - (f12 * f21)
            if D >= 0.0:
                Dmax = min(p * (1.0-q), q * (1.0-p))
            else:
                Dmax = min(p*q,(1-p)*(1-q))
            Dprime = D / Dmax
            rsquared = (D ** 2) / (p * (1-p) * q * (1-q))
            result["alphaDprime"] = round(Dprime,3)
            result["alpharsquared"] = round(rsquared,4)
            result["alphaf11"] = f11
            result["alphaf12"] = f12
            result["alphaf21"] = f21
            result["alphaf22"] = f22
            if min(f11,f12,f21,f22) < -0.0000001:
                result["alphaposs"] = 0
            else:
                result["alpha-e1111"] = n * f11**2
                result["alpha-e1112"] = 2 * n * f11 * f12
                result["alpha-e1122"] = n * f12**2
                result["alpha-e1211"] = 2 * n * f11 * f21
                result["alpha-e1212"] = 2 * n * f12 * f21 + 2 * n * f11 * f22 
                result["alpha-e1222"] = 2 * n * f12 * f22
                result["alpha-e2211"] = n * f21**2
                result["alpha-e2212"] = 2 * n * f21 * f22
                result["alpha-e2222"] = n * f22**2
                if result["alpha-e1111"] > 0.0:
                    chisq1111 = ((n1111 - result["alpha-e1111"])**2/result["alpha-e1111"])
                else: chisq1111 = 0.0
                if result["alpha-e1112"] > 0.0:
                    chisq1112 = ((n1112 - result["alpha-e1112"])**2/result["alpha-e1112"])
                else: chisq1112 = 0.0
                if result["alpha-e1122"] > 0.0:
                    chisq1122 = ((n1122 - result["alpha-e1122"])**2/result["alpha-e1122"])
                else: chisq1122 = 0.0
                if result["alpha-e1211"] > 0.0:
                    chisq1211 = ((n1211 - result["alpha-e1211"])**2/result["alpha-e1211"])
                else: chisq1211 = 0.0
                if result["alpha-e1212"] > 0.0:
                    chisq1212 = ((n1212 - result["alpha-e1212"])**2/result["alpha-e1212"])
                else: chisq1212 = 0.0
                if result["alpha-e1222"] > 0.0:
                    chisq1222 = ((n1222 - result["alpha-e1222"])**2/result["alpha-e1222"])
                else: chisq1222 = 0.0
                if result["alpha-e2211"] > 0.0:
                    chisq2211 = ((n2211 - result["alpha-e2211"])**2/result["alpha-e2211"])
                else: chisq2211 = 0.0
                if result["alpha-e2212"] > 0.0:
                    chisq2212 = ((n2212 - result["alpha-e2212"])**2/result["alpha-e2212"])
                else: chisq2212 = 0.0
                if result["alpha-e2222"] > 0.0:
                    chisq2222 = ((n2222 - result["alpha-e2222"])**2/result["alpha-e2222"])
                else: chisq2222 = 0.0
                result["alpha-chisq"] = (chisq1111+chisq1112+chisq1122+chisq1211+chisq1212+chisq1222+chisq2211+chisq2212+chisq2222)
        else:
            result["alphaposs"] = 0
        result["betaposs"] = 0
        result["gammaposs"] = 0

    elif yN2 == h2 and result["realnoproblem"] == 2: # Deactivated by making 2
        # option 2
        delta = math.pow((yN/2.0*a),(1.0/3.0))
        result["inputdata"] = str(inputarray)
        result["alpha"] = xN + delta
        result["beta"] = xN + delta
        result["gamma"] = xN - 2.0*delta
        result["p"] = str(p)
        result["q"] = str(q)
        result["pq"] = str(p * q)
        if result["alpha"] >= minhap - 0.00001 and result["alpha"] <= maxhap + 0.00001:
            result["alphaposs"] = 1
            f11 = result["alpha"]
            f12 = p - f11
            f21 = q - f11
            f22 = 1 - (f11 + f12 + f21)
            D = (f11 * f22) - (f12 * f21)
            if D >= 0.0:
                Dmax = min(p * (1.0-q), q * (1.0-p))
            else:
                Dmax = min(p*q,(1-p)*(1-q))
            Dprime = D / Dmax
            rsquared = (D ** 2) / (p * (1-p) * q * (1-q))
            result["alphaDprime"] = round(Dprime,3)
            result["alpharsquared"] = round(rsquared,4)
            result["alphaf11"] = f11
            result["alphaf12"] = f12
            result["alphaf21"] = f21
            result["alphaf22"] = f22
            if min(f11,f12,f21,f22) < -0.0000001:
                result["alphaposs"] = 0
            else:
                result["alpha-e1111"] = n * f11**2
                result["alpha-e1112"] = 2 * n * f11 * f12
                result["alpha-e1122"] = n * f12**2
                result["alpha-e1211"] = 2 * n * f11 * f21
                result["alpha-e1212"] = 2 * n * f12 * f21 + 2 * n * f11 * f22 
                result["alpha-e1222"] = 2 * n * f12 * f22
                result["alpha-e2211"] = n * f21**2
                result["alpha-e2212"] = 2 * n * f21 * f22
                result["alpha-e2222"] = n * f22**2
                if result["alpha-e1111"] > 0.0:
                    chisq1111 = ((n1111 - result["alpha-e1111"])**2/result["alpha-e1111"])
                else: chisq1111 = 0.0
                if result["alpha-e1112"] > 0.0:
                    chisq1112 = ((n1112 - result["alpha-e1112"])**2/result["alpha-e1112"])
                else: chisq1112 = 0.0
                if result["alpha-e1122"] > 0.0:
                    chisq1122 = ((n1122 - result["alpha-e1122"])**2/result["alpha-e1122"])
                else: chisq1122 = 0.0
                if result["alpha-e1211"] > 0.0:
                    chisq1211 = ((n1211 - result["alpha-e1211"])**2/result["alpha-e1211"])
                else: chisq1211 = 0.0
                if result["alpha-e1212"] > 0.0:
                    chisq1212 = ((n1212 - result["alpha-e1212"])**2/result["alpha-e1212"])
                else: chisq1212 = 0.0
                if result["alpha-e1222"] > 0.0:
                    chisq1222 = ((n1222 - result["alpha-e1222"])**2/result["alpha-e1222"])
                else: chisq1222 = 0.0
                if result["alpha-e2211"] > 0.0:
                    chisq2211 = ((n2211 - result["alpha-e2211"])**2/result["alpha-e2211"])
                else: chisq2211 = 0.0
                if result["alpha-e2212"] > 0.0:
                    chisq2212 = ((n2212 - result["alpha-e2212"])**2/result["alpha-e2212"])
                else: chisq2212 = 0.0
                if result["alpha-e2222"] > 0.0:
                    chisq2222 = ((n2222 - result["alpha-e2222"])**2/result["alpha-e2222"])
                else: chisq2222 = 0.0
                result["alpha-chisq"] = (chisq1111+chisq1112+chisq1122+chisq1211+chisq1212+chisq1222+chisq2211+chisq2212+chisq2222)
        else:
            result["alphaposs"] = 0
        if result["beta"] >= minhap - 0.00001 and result["beta"] <= maxhap + 0.00001:
            result["betaposs"] = 1
            f11 = result["beta"]
            f12 = p - f11
            f21 = q - f11
            f22 = 1 - (f11 + f12 + f21)
            D = (f11 * f22) - (f12 * f21)
            if D >= 0.0:
                Dmax = min(p * (1.0-q), q * (1.0-p))
            else:
                Dmax = min(p*q,(1-p)*(1-q))
            Dprime = D / Dmax
            rsquared = (D ** 2) / (p * (1-p) * q * (1-q))
            result["betaDprime"] = round(Dprime,3)
            result["betarsquared"] = round(rsquared,4)
            result["betaf11"] = f11
            result["betaf12"] = f12
            result["betaf21"] = f21
            result["betaf22"] = f22
            if min(f11,f12,f21,f22) < -0.0000001:
                result["betaposs"] = 0
            else:
                result["beta-e1111"] = n * f11**2
                result["beta-e1112"] = 2 * n * f11 * f12
                result["beta-e1122"] = n * f12**2
                result["beta-e1211"] = 2 * n * f11 * f21
                result["beta-e1212"] = 2 * n * f12 * f21 + 2 * n * f11 * f22 
                result["beta-e1222"] = 2 * n * f12 * f22
                result["beta-e2211"] = n * f21**2
                result["beta-e2212"] = 2 * n * f21 * f22
                result["beta-e2222"] = n * f22**2
                if result["beta-e1111"] > 0.0:
                    chisq1111 = ((n1111 - result["beta-e1111"])**2/result["beta-e1111"])
                else: chisq1111 = 0.0
                if result["beta-e1112"] > 0.0:
                    chisq1112 = ((n1112 - result["beta-e1112"])**2/result["beta-e1112"])
                else: chisq1112 = 0.0
                if result["beta-e1122"] > 0.0:
                    chisq1122 = ((n1122 - result["beta-e1122"])**2/result["beta-e1122"])
                else: chisq1122 = 0.0
                if result["beta-e1211"] > 0.0:
                    chisq1211 = ((n1211 - result["beta-e1211"])**2/result["beta-e1211"])
                else: chisq1211 = 0.0
                if result["beta-e1212"] > 0.0:
                    chisq1212 = ((n1212 - result["beta-e1212"])**2/result["beta-e1212"])
                else: chisq1212 = 0.0
                if result["beta-e1222"] > 0.0:
                    chisq1222 = ((n1222 - result["beta-e1222"])**2/result["beta-e1222"])
                else: chisq1222 = 0.0
                if result["beta-e2211"] > 0.0:
                    chisq2211 = ((n2211 - result["beta-e2211"])**2/result["beta-e2211"])
                else: chisq2211 = 0.0
                if result["beta-e2212"] > 0.0:
                    chisq2212 = ((n2212 - result["beta-e2212"])**2/result["beta-e2212"])
                else: chisq2212 = 0.0
                if result["beta-e2222"] > 0.0:
                    chisq2222 = ((n2222 - result["beta-e2222"])**2/result["beta-e2222"])
                else: chisq2222 = 0.0
                result["beta-chisq"] = (chisq1111+chisq1112+chisq1122+chisq1211+chisq1212+chisq1222+chisq2211+chisq2212+chisq2222)
        else:
            result["betaposs"] = 0
        if result["gamma"] >= minhap - 0.00001 and result["gamma"] <= maxhap + 0.00001:
            result["gammaposs"] = 1
            f11 = result["gamma"]
            f12 = p - f11
            f21 = q - f11
            f22 = 1 - (f11 + f12 + f21)
            D = (f11 * f22) - (f12 * f21)
            if D >= 0.0:
                Dmax = min(p * (1.0-q), q * (1.0-p))
            else:
                Dmax = min(p*q,(1-p)*(1-q))
            Dprime = D / Dmax
            rsquared = (D ** 2) / (p * (1-p) * q * (1-q))
            result["gammaDprime"] = round(Dprime,3)
            result["gammarsquared"] = round(rsquared,4)
            result["gammaf11"] = f11
            result["gammaf12"] = f12
            result["gammaf21"] = f21
            result["gammaf22"] = f22
            if min(f11,f12,f21,f22) < -0.0000001:
                result["gammaposs"] = 0
            else:
                result["gamma-e1111"] = n * f11**2
                result["gamma-e1112"] = 2 * n * f11 * f12
                result["gamma-e1122"] = n * f12**2
                result["gamma-e1211"] = 2 * n * f11 * f21
                result["gamma-e1212"] = 2 * n * f12 * f21 + 2 * n * f11 * f22 
                result["gamma-e1222"] = 2 * n * f12 * f22
                result["gamma-e2211"] = n * f21**2
                result["gamma-e2212"] = 2 * n * f21 * f22
                result["gamma-e2222"] = n * f22**2
                if result["gamma-e1111"] > 0.0:
                    chisq1111 = ((n1111 - result["gamma-e1111"])**2/result["gamma-e1111"])
                else: chisq1111 = 0.0
                if result["gamma-e1112"] > 0.0:
                    chisq1112 = ((n1112 - result["gamma-e1112"])**2/result["gamma-e1112"])
                else: chisq1112 = 0.0
                if result["gamma-e1122"] > 0.0:
                    chisq1122 = ((n1122 - result["gamma-e1122"])**2/result["gamma-e1122"])
                else: chisq1122 = 0.0
                if result["gamma-e1211"] > 0.0:
                    chisq1211 = ((n1211 - result["gamma-e1211"])**2/result["gamma-e1211"])
                else: chisq1211 = 0.0
                if result["gamma-e1212"] > 0.0:
                    chisq1212 = ((n1212 - result["gamma-e1212"])**2/result["gamma-e1212"])
                else: chisq1212 = 0.0
                if result["gamma-e1222"] > 0.0:
                    chisq1222 = ((n1222 - result["gamma-e1222"])**2/result["gamma-e1222"])
                else: chisq1222 = 0.0
                if result["gamma-e2211"] > 0.0:
                    chisq2211 = ((n2211 - result["gamma-e2211"])**2/result["gamma-e2211"])
                else: chisq2211 = 0.0
                if result["gamma-e2212"] > 0.0:
                    chisq2212 = ((n2212 - result["gamma-e2212"])**2/result["gamma-e2212"])
                else: chisq2212 = 0.0
                if result["gamma-e2222"] > 0.0:
                    chisq2222 = ((n2222 - result["gamma-e2222"])**2/result["gamma-e2222"])
                else: chisq2222 = 0.0
                result["gamma-chisq"] = (chisq1111+chisq1112+chisq1122+chisq1211+chisq1212+chisq1222+chisq2211+chisq2212+chisq2222)
        else:
            result["gammaposs"] = 0
        
    else:
        #option 3
        h = math.pow(h2, 0.5)
        theta = ((math.acos(-yN/h))/3.0)
        #delta = math.pow((yN/2.0*a),(1.0/3.0)) # is it correct to reuse this?
        delta = math.pow(d2,0.5)
        result["inputdata"] = inputarray
        result["alpha"] = xN + 2.0 * delta * math.cos(theta)
        result["beta"] = xN + 2.0 * delta * math.cos(2.0 * math.pi/3.0 + theta)
        result["gamma"] = xN + 2.0 * delta * math.cos(4.0 * math.pi/3.0 + theta)
        result["p"] = p
        result["q"] = q
        result["pq"] = p * q
        if result["alpha"] >= minhap - 0.00001 and result["alpha"] <= maxhap + 0.00001:
            result["alphaposs"] = 1
            f11 = result["alpha"]
            f12 = p - f11
            f21 = q - f11
            f22 = 1 - (f11 + f12 + f21)
            D = (f11 * f22) - (f12 * f21)
            if D >= 0.0:
                Dmax = min(p * (1.0-q), q * (1.0-p))
            else:
                Dmax = min(p*q,(1-p)*(1-q))
            Dprime = D / Dmax
            rsquared = (D ** 2) / (p * (1-p) * q * (1-q))
            result["alphaDprime"] = round(Dprime,3)
            result["alpharsquared"] = round(rsquared,4)
            result["alphaf11"] = f11
            result["alphaf12"] = f12
            result["alphaf21"] = f21
            result["alphaf22"] = f22
            if min(f11,f12,f21,f22) < -0.0000001:
                result["alphaposs"] = 0
            else:
                result["alpha-e1111"] = n * f11**2
                result["alpha-e1112"] = 2 * n * f11 * f12
                result["alpha-e1122"] = n * f12**2
                result["alpha-e1211"] = 2 * n * f11 * f21
                result["alpha-e1212"] = 2 * n * f12 * f21 + 2 * n * f11 * f22 
                result["alpha-e1222"] = 2 * n * f12 * f22
                result["alpha-e2211"] = n * f21**2
                result["alpha-e2212"] = 2 * n * f21 * f22
                result["alpha-e2222"] = n * f22**2
                if result["alpha-e1111"] > 0.0:
                    chisq1111 = ((n1111 - result["alpha-e1111"])**2/result["alpha-e1111"])
                else: chisq1111 = 0.0
                if result["alpha-e1112"] > 0.0:
                    chisq1112 = ((n1112 - result["alpha-e1112"])**2/result["alpha-e1112"])
                else: chisq1112 = 0.0
                if result["alpha-e1122"] > 0.0:
                    chisq1122 = ((n1122 - result["alpha-e1122"])**2/result["alpha-e1122"])
                else: chisq1122 = 0.0
                if result["alpha-e1211"] > 0.0:
                    chisq1211 = ((n1211 - result["alpha-e1211"])**2/result["alpha-e1211"])
                else: chisq1211 = 0.0
                if result["alpha-e1212"] > 0.0:
                    chisq1212 = ((n1212 - result["alpha-e1212"])**2/result["alpha-e1212"])
                else: chisq1212 = 0.0
                if result["alpha-e1222"] > 0.0:
                    chisq1222 = ((n1222 - result["alpha-e1222"])**2/result["alpha-e1222"])
                else: chisq1222 = 0.0
                if result["alpha-e2211"] > 0.0:
                    chisq2211 = ((n2211 - result["alpha-e2211"])**2/result["alpha-e2211"])
                else: chisq2211 = 0.0
                if result["alpha-e2212"] > 0.0:
                    chisq2212 = ((n2212 - result["alpha-e2212"])**2/result["alpha-e2212"])
                else: chisq2212 = 0.0
                if result["alpha-e2222"] > 0.0:
                    chisq2222 = ((n2222 - result["alpha-e2222"])**2/result["alpha-e2222"])
                else: chisq2222 = 0.0
                result["alpha-chisq"] = (chisq1111+chisq1112+chisq1122+chisq1211+chisq1212+chisq1222+chisq2211+chisq2212+chisq2222)
        else: result["alphaposs"] = 0
        if result["beta"] >= minhap - 0.00001 and result["beta"] <= maxhap + 0.00001:
            result["betaposs"] = 1
            f11 = result["beta"]
            f12 = p - f11
            f21 = q - f11
            f22 = 1 - (f11 + f12 + f21)
            D = (f11 * f22) - (f12 * f21)
            if D >= 0.0:
                Dmax = min(p * (1.0-q), q * (1.0-p))
            else:
                Dmax = min(p*q,(1-p)*(1-q))
            Dprime = D / Dmax
            rsquared = (D ** 2) / (p * (1-p) * q * (1-q))
            result["betaDprime"] = round(Dprime,3)
            result["betarsquared"] = round(rsquared,4)
            result["betaf11"] = f11
            result["betaf12"] = f12
            result["betaf21"] = f21
            result["betaf22"] = f22
            if min(f11,f12,f21,f22) < -0.0000001:
                result["betaposs"] = 0
            else:
                result["beta-e1111"] = n * f11**2
                result["beta-e1112"] = 2 * n * f11 * f12
                result["beta-e1122"] = n * f12**2
                result["beta-e1211"] = 2 * n * f11 * f21
                result["beta-e1212"] = 2 * n * f12 * f21 + 2 * n * f11 * f22
                result["beta-e1222"] = 2 * n * f12 * f22
                result["beta-e2211"] = n * f21**2
                result["beta-e2212"] = 2 * n * f21 * f22
                result["beta-e2222"] = n * f22**2
                if result["beta-e1111"] > 0.0:
                    chisq1111 = ((n1111 - result["beta-e1111"])**2/result["beta-e1111"])
                else: chisq1111 = 0.0
                if result["beta-e1112"] > 0.0:
                    chisq1112 = ((n1112 - result["beta-e1112"])**2/result["beta-e1112"])
                else: chisq1112 = 0.0
                if result["beta-e1122"] > 0.0:
                    chisq1122 = ((n1122 - result["beta-e1122"])**2/result["beta-e1122"])
                else: chisq1122 = 0.0
                if result["beta-e1211"] > 0.0:
                    chisq1211 = ((n1211 - result["beta-e1211"])**2/result["beta-e1211"])
                else: chisq1211 = 0.0
                if result["beta-e1212"] > 0.0:
                    chisq1212 = ((n1212 - result["beta-e1212"])**2/result["beta-e1212"])
                else: chisq1212 = 0.0
                if result["beta-e1222"] > 0.0:
                    chisq1222 = ((n1222 - result["beta-e1222"])**2/result["beta-e1222"])
                else: chisq1222 = 0.0
                if result["beta-e2211"] > 0.0:
                    chisq2211 = ((n2211 - result["beta-e2211"])**2/result["beta-e2211"])
                else: chisq2211 = 0.0
                if result["beta-e2212"] > 0.0:
                    chisq2212 = ((n2212 - result["beta-e2212"])**2/result["beta-e2212"])
                else: chisq2212 = 0.0
                if result["beta-e2222"] > 0.0:
                    chisq2222 = ((n2222 - result["beta-e2222"])**2/result["beta-e2222"])
                else: chisq2222 = 0.0
                result["beta-chisq"] = (chisq1111+chisq1112+chisq1122+chisq1211+chisq1212+chisq1222+chisq2211+chisq2212+chisq2222)
        else: result["betaposs"] = 0
        if result["gamma"] >= minhap - 0.00001 and result["gamma"] <= maxhap + 0.00001:
            result["gammaposs"] = 1
            f11 = result["gamma"]
            f12 = p - f11
            f21 = q - f11
            f22 = 1 - (f11 + f12 + f21)
            D = (f11 * f22) - (f12 * f21)
            if D >= 0.0:
                Dmax = min(p * (1.0-q), q * (1.0-p))
            else:
                Dmax = min(p*q,(1-p)*(1-q))
            Dprime = D / Dmax
            rsquared = (D ** 2) / (p * (1-p) * q * (1-q))
            result["gammaDprime"] = round(Dprime,3)
            result["gammarsquared"] = round(rsquared,4)
            result["gammaf11"] = f11
            result["gammaf12"] = f12
            result["gammaf21"] = f21
            result["gammaf22"] = f22
            if min(f11,f12,f21,f22) < -0.0000001:
                result["gammaposs"] = 0
            else:
                result["gamma-e1111"] = n * f11**2 
                result["gamma-e1112"] = 2 * n * f11 * f12
                result["gamma-e1122"] = n * f12**2
                result["gamma-e1211"] = 2 * n * f11 * f21
                result["gamma-e1212"] = 2 * n * f12 * f21 + 2 * n * f11 * f22 
                result["gamma-e1222"] = 2 * n * f12 * f22
                result["gamma-e2211"] = n * f21**2
                result["gamma-e2212"] = 2 * n * f21 * f22
                result["gamma-e2222"] = n * f22**2
                if result["gamma-e1111"] > 0.0:
                    chisq1111 = ((n1111 - result["gamma-e1111"])**2/result["gamma-e1111"])
                else: chisq1111 = 0.0
                if result["gamma-e1112"] > 0.0:
                    chisq1112 = ((n1112 - result["gamma-e1112"])**2/result["gamma-e1112"])
                else: chisq1112 = 0.0
                if result["gamma-e1122"] > 0.0:
                    chisq1122 = ((n1122 - result["gamma-e1122"])**2/result["gamma-e1122"])
                else: chisq1122 = 0.0
                if result["gamma-e1211"] > 0.0:
                    chisq1211 = ((n1211 - result["gamma-e1211"])**2/result["gamma-e1211"])
                else: chisq1211 = 0.0
                if result["gamma-e1212"] > 0.0:
                    chisq1212 = ((n1212 - result["gamma-e1212"])**2/result["gamma-e1212"])
                else: chisq1212 = 0.0
                if result["gamma-e1222"] > 0.0:
                    chisq1222 = ((n1222 - result["gamma-e1222"])**2/result["gamma-e1222"])
                else: chisq1222 = 0.0
                if result["gamma-e2211"] > 0.0:
                    chisq2211 = ((n2211 - result["gamma-e2211"])**2/result["gamma-e2211"])
                else: chisq2211 = 0.0
                if result["gamma-e2212"] > 0.0:
                    chisq2212 = ((n2212 - result["gamma-e2212"])**2/result["gamma-e2212"])
                else: chisq2212 = 0.0
                if result["gamma-e2222"] > 0.0:
                    chisq2222 = ((n2222 - result["gamma-e2222"])**2/result["gamma-e2222"])
                else: chisq2222 = 0.0
                result["gamma-chisq"] = (chisq1111+chisq1112+chisq1122+chisq1211+chisq1212+chisq1222+chisq2211+chisq2212+chisq2222)
        else: result["gammaposs"] = 0
        
    #else: result["Error"] = "No answer"
    # Return the result as a dictionary array
    return result

n1111 = float(form["n1111"].value)
n1112 = float(form["n1112"].value)
n1122 = float(form["n1122"].value)
n1211 = float(form["n1211"].value)
n1212 = float(form["n1212"].value)
n1222 = float(form["n1222"].value)
n2211 = float(form["n2211"].value)
n2212 = float(form["n2212"].value)
n2222 = float(form["n2222"].value)
result = cubex([n1111,n1112,n1122,n1211,n1212,n1222,n2211,n2212,n2222])

print "Content-type: text/html\n\n";
print """<html><head><style>p{margin:8px;}</style></head><body BGCOLOR="#ebebeb"><div style="margin: 0 auto; width: 802px; background-color: rgb(255,255,255); text-align: centre; border: 1px solid rgb(187,187,187);"> """
if result["realnoproblem"] == 2: # Deactivated by making 2
    print """<div style="width: 800px; text-align: center;" id="layer1">
<img style="width: 585px; height: 75px;" alt=""
 src="http://www.oege.org/software/cubex/cubex-title.jpg"></div>
<H1 style="border-top: 1px dashed rgb(187, 187, 187); background: rgb(255, 255, 255) none repeat scroll 0% 50%; -moz-background-clip: initial; -moz-background-origin: initial; -moz-background-inline-policy: initial; text-align: center; position: relative; width: 800px; z-index: 1;">Err 1: There is a problem with real number calculation for his dataset</H1>Please contact cubex@genes.org.uk with details of the problem or any queries."""
else:
    print """<div style="width: 800px; text-align: center;" id="layer1">
<img style="width: 585px; height: 75px;" alt=""
 src="../cubex-title.jpg"></div>
<h3 style="border-top: 1px dashed rgb(187, 187, 187); background: rgb(255, 255, 255) none repeat scroll 0% 50%; -moz-background-clip: initial; -moz-background-origin: initial; -moz-background-inline-policy: initial; text-align: center; position: relative; width: 800px; z-index: 1;" id="layer1">Results </H3>"""
    print """<p>For an explanation of the analysis and results please see <a href="#note">notes</a> below.</p>"""
    if result["hwchisnp1"] >= 3.84 and result["hwchisnp2"] >= 3.84:
        print """<p><div style="background:#FFFF55;z-index: 1;margin: 8px;" id="layer1"><font style="color:#FF0000"><b>Both SNPs are significantly out of Hardy-Weinberg Equilibrium</b> - this is inconsistent with the assumptions of the model. Please use <a href="http://www.oege.org/software/hardy-weinberg.shtml">this calculator</a> to check your data.</font></div></p>"""
    elif result["hwchisnp1"] >= 3.84 and result["hwchisnp2"] < 3.84:
        print """<p><div style="background:#FFFF55;z-index: 1;margin: 8px;" id="layer1"><font style="color:#FF0000"><b>SNP 1 is significantly out of Hardy-Weinberg Equilibrium</b> - this is inconsistent with the assumptions of the model. Please use <a href="http://www.oege.org/software/hardy-weinberg.shtml">this calculator</a> to check your data.</font></div></p>"""
    elif result["hwchisnp2"] >= 3.84 and result["hwchisnp1"] < 3.84:
        print """<p><div style="background:#FFFF55;z-index: 1;margin: 8px;" id="layer1"><font style="color:#FF0000"><b>SNP 2 is significantly out of Hardy-Weinberg Equilibrium</b> - this is inconsistent with the assumptions of the model. Please use <a href="http://www.oege.org/software/hardy-weinberg.shtml">this calculator</a> to check your data.</font></div></p>"""

    print """<p><b>Number of biologically possible solutions: </b>"""
    print str(result["alphaposs"] + result["betaposs"] + result["gammaposs"])
    print """. """
    if result["alphaposs"] == 1 and result["betaposs"] == 1 and result["gammaposs"] == 1:
        if result["alpha-chisq"] == min(result["alpha-chisq"],result["beta-chisq"],result["gamma-chisq"]):
            print """<br><b>&#945;</b> is the most likely solution<br>"""
        if result["beta-chisq"] == min(result["alpha-chisq"],result["beta-chisq"],result["gamma-chisq"]):
            print """<br><b>&#946;</b> is the most likely solution<br>"""
        if result["gamma-chisq"] == min(result["alpha-chisq"],result["beta-chisq"],result["gamma-chisq"]):
            print """<br><b>&#947;</b> is the most likely solution<br>"""
    if result["alphaposs"] == 1 and result["betaposs"] == 1 and result["gammaposs"] == 0:
        if result["alpha-chisq"] == min(result["alpha-chisq"],result["beta-chisq"]):
            print """<br><b>&#945;</b> is the most likely solution<br>"""
        if result["beta-chisq"] == min(result["alpha-chisq"],result["beta-chisq"]):
            print """<br><b>&#946;</b> is the most likely solution<br>"""
    if result["alphaposs"] == 1 and result["betaposs"] == 0 and result["gammaposs"] == 1:
        if result["alpha-chisq"] == min(result["alpha-chisq"],result["gamma-chisq"]):
            print """<br><b>&#945;</b> is the most likely solution<br>"""
        if result["gamma-chisq"] == min(result["alpha-chisq"],result["gamma-chisq"]):
            print """<br><b>&#947;</b> is the most likely solution<br>"""
    if result["alphaposs"] == 0 and result["betaposs"] == 1 and result["gammaposs"] == 1:
        if result["beta-chisq"] == min(result["beta-chisq"],result["gamma-chisq"]):
            print """<br><b>&#946;</b> is the most likely solution<br>"""
        if result["gamma-chisq"] == min(result["beta-chisq"],result["gamma-chisq"]):
            print """<br><b>&#947;</b> is the most likely solution<br>"""
    print """See <a href="#chisq"><i>X</i><sup>2</sup> table</a> below 3x3</p><br>"""
    print """<table border="1" cellspacing="0" bordercolordark="#000000" align="center" bordercolorlight="#FFFFFF"><tr><td></td><td colspan=4 align="center"><b>Haplotype frequencies</b></td><td width=10></td><td colspan=3 align="center"><b>LD statistics</b></td></tr>"""
    print """<tr><td width=65 align="center">Solution</td><td width=65 align="center">&#402;<sub>11</sub></td><td width=65 align="center">&#402;<sub>12</sub></td><td width=65 align="center">&#402;<sub>21</sub></td>"""
    print """<td width=65 align="center">&#402;<sub>22</sub></td><td></td><td width=65 align="center">D'</sup></td><td width=65 align="center">r<sup>2</sup></td><td width=65 align="center"><i>X</i><sup>2</sup></td></tr>"""
    if result["alphaposs"] == 1:
        print """<tr><td align="center"><b>&#945;</b></td><td align="center"><div style="background:#FFCCCC;z-index: 1;color:#FF0000" id="layer1">"""
        print round(result["alphaf11"],4)
        print """</div></td><td align="center"><div style="background:#FFCCCC;z-index: 1;color:#FF0000" id="layer1">"""
        print round(result["alphaf12"],4)
        print """</div></td><td align="center"><div style="background:#FFCCCC;z-index: 1;color:#FF0000" id="layer1">"""
        print round(result["alphaf21"],4)
        print """</div></td><td align="center"><div style="background:#FFCCCC;z-index: 1;color:#FF0000" id="layer1">"""
        print round(result["alphaf22"],4)
        print """</div></td><td></td><td align="center"><div style="background:#FFCCCC;z-index: 1;color:#FF0000" id="layer1">"""
        print result["alphaDprime"]
        print """</div></td><td align="center"><div style="background:#FFCCCC;z-index: 1;color:#FF0000" id="layer1">"""
        print result["alpharsquared"]
        print """</div></td><td align="center"><div style="background:#FFCCCC;z-index: 1;color:#FF0000" id="layer1">"""
        print round(float(result["alpharsquared"]) * float(result["n"]),2)
        print """</div></td></tr>"""
    if result["betaposs"] == 1:
        print """<tr><td align="center"><b>&#946;</b></td><td align="center"><div style="background:#CCFFCC;z-index: 1;color:#00FF00" id="layer1">"""
        print round(result["betaf11"],4)
        print """</div></td><td align="center"><div style="background:#CCFFCC;z-index: 1;color:#00FF00" id="layer1">"""
        print round(result["betaf12"],4)
        print """</div></td><td align="center"><div style="background:#CCFFCC;z-index: 1;color:#00FF00" id="layer1">"""
        print round(result["betaf21"],4)
        print """</div></td><td align="center"><div style="background:#CCFFCC;z-index: 1;color:#00FF00" id="layer1">"""
        print round(result["betaf22"],4)
        print """</div></td><td></td><td align="center"><div style="background:#CCFFCC;z-index: 1;color:#00FF00" id="layer1">"""
        print result["betaDprime"]
        print """</div></td><td align="center"><div style="background:#CCFFCC;z-index: 1;color:#00FF00" id="layer1">"""
        print result["betarsquared"]
        print """</div></td><td align="center"><div style="background:#CCFFCC;z-index: 1;color:#00FF00" id="layer1">"""
        print round(float(result["betarsquared"]) * float(result["n"]),2)
        print """</div></td></tr>"""
    if result["gammaposs"] == 1:
        print """<tr><td align="center"><b>&#947;</b></td><td align="center"><div style="background:#CCCCFF;z-index: 1;color:#0000FF" id="layer1">"""
        print round(result["gammaf11"],4)
        print """</div></td><td align="center"><div style="background:#CCCCFF;z-index: 1;color:#0000FF" id="layer1">"""
        print round(result["gammaf12"],4)
        print """</div></td><td align="center"><div style="background:#CCCCFF;z-index: 1;color:#0000FF" id="layer1">"""
        print round(result["gammaf21"],4)
        print """</div></td><td align="center"><div style="background:#CCCCFF;z-index: 1;color:#0000FF" id="layer1">"""
        print round(result["gammaf22"],4)
        print """</div></td><td></td><td align="center"><div style="background:#CCCCFF;z-index: 1;color:#0000FF" id="layer1">"""
        print result["gammaDprime"]
        print """</div></td><td align="center"><div style="background:#CCCCFF;z-index: 1;color:#0000FF" id="layer1">"""
        print result["gammarsquared"]
        print """</div></td><td align="center"><div style="background:#CCCCFF;z-index: 1;color:#0000FF" id="layer1">"""
        print round(float(result["gammarsquared"]) * float(result["n"]),2)
        print """</div></td></tr>"""
    print """</table>"""
    print """<br>"""

    print """<H3 style="border-top: 1px dashed rgb(187, 187, 187); background: rgb(255, 255, 255) none repeat scroll 0% 50%; -moz-background-clip: initial; -moz-background-origin: initial; -moz-background-inline-policy: initial; text-align: center; position: relative; width: 800px; z-index: 1;">3x3 table of observed and expected diplotype numbers</H3>"""
    print """<p>Black numbers on white are original data entered<br>Coloured numbers on coloured background represent """
    print """the solutions from the table above.</p> """
    print """<table border="1" align="center" width="400" bordercolorlight="#C0C0C0" cellspacing="0" bordercolordark="#000000"><tr><td align="center" width="115" bordercolorlight="#FFFFFF" bordercolordark="#FFFFFF" colspan="2" rowspan="2"></td><td colspan="3" align="center"><b>SNP 2</b></td></tr><tr><td align="center">11</td><td align="center">12</td><td align="center">22</td></tr><tr><td rowspan="3" align="center" width="80"><b>SNP 1</b></td><td align="center" width="51">11</td><td align="center">"""
    print int(n1111)
    if result["alphaposs"] == 1:
        print """<div style="background:#FFCCCC;z-index: 1;color:#FF0000" id="layer1">"""
        print round(result["alpha-e1111"],1)
        print """</div>"""
    if result["betaposs"] == 1:
        print """<div style="background:#CCFFCC;z-index: 1;color:#00FF00" id="layer1">"""
        print round(result["beta-e1111"],1)
        print """</div>"""
    if result["gammaposs"] == 1:
        print """<div style="background:#CCCCFF;z-index: 1;color:#0000FF" id="layer1">"""
        print round(result["gamma-e1111"],1)
        print """</div>"""
    print """</td>	<td align="center">"""
    print int(n1112)
    if result["alphaposs"] == 1:
        print """<div style="background:#FFCCCC;z-index: 1;color:#FF0000" id="layer1">"""
        print round(result["alpha-e1112"],1)
        print """</div>"""
    if result["betaposs"] == 1:
        print """<div style="background:#CCFFCC;z-index: 1;color:#00FF00" id="layer1">"""
        print round(result["beta-e1112"],1)
        print """</div>"""
    if result["gammaposs"] == 1:
        print """<div style="background:#CCCCFF;z-index: 1;color:#0000FF" id="layer1">"""
        print round(result["gamma-e1112"],1)
        print """</div>"""
    print """</td>	<td align="center">"""
    print int(n1122)
    if result["alphaposs"] == 1:
        print """<div style="background:#FFCCCC;z-index: 1;color:#FF0000" id="layer1">"""
        print round(result["alpha-e1122"],1)
        print """</div>"""
    if result["betaposs"] == 1:
        print """<div style="background:#CCFFCC;z-index: 1;color:#00FF00" id="layer1">"""
        print round(result["beta-e1122"],1)
        print """</div>"""
    if result["gammaposs"] == 1:
        print """<div style="background:#CCCCFF;z-index: 1;color:#0000FF" id="layer1">"""
        print round(result["gamma-e1122"],1)
        print """</div>"""
    print """</td></tr><tr><td align="center" width="51">12</td><td align="center">"""
    print int(n1211)
    if result["alphaposs"] == 1:
        print """<div style="background:#FFCCCC;z-index: 1;color:#FF0000" id="layer1">"""
        print round(result["alpha-e1211"],1)
        print """</div>"""
    if result["betaposs"] == 1:
        print """<div style="background:#CCFFCC;z-index: 1;color:#00FF00" id="layer1">"""
        print round(result["beta-e1211"],1)
        print """</div>"""
    if result["gammaposs"] == 1:
        print """<div style="background:#CCCCFF;z-index: 1;color:#0000FF" id="layer1">"""
        print round(result["gamma-e1211"],1)
        print """</div>"""
    print """</td><td align="center">"""
    print int(n1212)
    if result["alphaposs"] == 1:
        print """<div style="background:#FFCCCC;z-index: 1;color:#FF0000" id="layer1">"""
        print round(result["alpha-e1212"],1)
        print """</div>"""
    if result["betaposs"] == 1:
        print """<div style="background:#CCFFCC;z-index: 1;color:#00FF00" id="layer1">"""
        print round(result["beta-e1212"],1)
        print """</div>"""
    if result["gammaposs"] == 1:
        print """<div style="background:#CCCCFF;z-index: 1;color:#0000FF" id="layer1">"""
        print round(result["gamma-e1212"],1)
        print """</div>"""
    print """</td><td align="center">"""
    print int(n1222)
    if result["alphaposs"] == 1:
        print """<div style="background:#FFCCCC;z-index: 1;color:#FF0000" id="layer1">"""
        print round(result["alpha-e1222"],1)
        print """</div>"""
    if result["betaposs"] == 1:
        print """<div style="background:#CCFFCC;z-index: 1;color:#00FF00" id="layer1">"""
        print round(result["beta-e1222"],1)
        print """</div>"""
    if result["gammaposs"] == 1:
        print """<div style="background:#CCCCFF;z-index: 1;color:#0000FF" id="layer1">"""
        print round(result["gamma-e1222"],1)
        print """</div>"""
    print """</td></tr><tr><td align="center" width="51">22</td><td align="center">"""
    print int(n2211)
    if result["alphaposs"] == 1:
        print """<div style="background:#FFCCCC;z-index: 1;color:#FF0000" id="layer1">"""
        print round(result["alpha-e2211"],1)
        print """</div>"""
    if result["betaposs"] == 1:
        print """<div style="background:#CCFFCC;z-index: 1;color:#00FF00" id="layer1">"""
        print round(result["beta-e2211"],1)
        print """</div>"""
    if result["gammaposs"] == 1:
        print """<div style="background:#CCCCFF;z-index: 1;color:#0000FF" id="layer1">"""
        print round(result["gamma-e2211"],1)
        print """</div>"""
    print """</td><td align="center">"""
    print int(n2212)
    if result["alphaposs"] == 1:
        print """<div style="background:#FFCCCC;z-index: 1;color:#FF0000" id="layer1">"""
        print round(result["alpha-e2212"],1)
        print """</div>"""
    if result["betaposs"] == 1:
        print """<div style="background:#CCFFCC;z-index: 1;color:#00FF00" id="layer1">"""
        print round(result["beta-e2212"],1)
        print """</div>"""
    if result["gammaposs"] == 1:
        print """<div style="background:#CCCCFF;z-index: 1;color:#0000FF" id="layer1">"""
        print round(result["gamma-e2212"],1)
        print """</div>"""
    print """</td><td align="center">"""
    print int(n2222)
    if result["alphaposs"] == 1:
        print """<div style="background:#FFCCCC;z-index: 1;color:#FF0000" id="layer1">"""
        print round(result["alpha-e2222"],1)
        print """</div>"""
    if result["betaposs"] == 1:
        print """<div style="background:#CCFFCC;z-index: 1;color:#00FF00" id="layer1">"""
        print round(result["beta-e2222"],1)
        print """</div>"""
    if result["gammaposs"] == 1:
        print """<div style="background:#CCCCFF;z-index: 1;color:#0000FF" id="layer1">"""
        print round(result["gamma-e2222"],1)
        print """</div>"""
    print """</td></tr></table>"""
    # Chi-square
    print """<br><table border="1" align="center" cellspacing="0" bordercolordark="#000000" bordercolorlight="#FFFFFF">"""
    print """<tr><td width=65 align="center"><a name="chisq">Solution</a></td>"""
    print """<td width=65 align="center"><i>&#935;</i><sup>2</sup> of 3x3</tr>"""
    if result["alphaposs"] == 1:
        print """<tr><td align="center"><b>&#945;</b></td><td align="center"><div style="background:#FFCCCC;z-index: 1;color:#FF0000" id="layer1">"""

        print round(result["alpha-chisq"],4)

        print """</div></td></tr>"""
    if result["betaposs"] == 1:
        print """<tr><td align="center"><b>&#946;</b></td><td align="center"><div style="background:#CCFFCC;z-index: 1;color:#00FF00" id="layer1">"""

        print round(result["beta-chisq"],4)

        print """</div></td></tr>"""
    if result["gammaposs"] == 1:
        print """<tr><td align="center"><b>&#947;</b></td><td align="center"><div style="background:#CCCCFF;z-index: 1;color:#0000FF" id="layer1">"""

        print round(result["gamma-chisq"],4)

        print """</div></td></tr>"""
    print """</table><br>"""
    print """<div style="background: #cccccc; width:784px; margin: 8px;">This is a <i>&#935;</i><sup>2</sup> of the """
    print """3x3 table. The higher the value, the less good the fit of the observed haplotypes to Hardy-Weinberg equilibrium. Please see footnote regarding degrees of freedom. """
    print """If there are two or more solutions, the lower values are more likely (although note the different degrees of freedom if there are empty cells). However, a significant value indicates genotype data out of Hardy-Weinberg equilibrium, a problem that should be addressed before interpreting these results.</div>"""
    # End chisquare

    print """<H3 style="border-top: 1px dashed rgb(187, 187, 187); background: rgb(255, 255, 255) none repeat scroll 0% 50%; -moz-background-clip: initial; -moz-background-origin: initial; -moz-background-inline-policy: initial; text-align: center; position: relative; width: 800px; z-index: 1;">Other statistics</H3>"""

    print """<p>Minimum <i>biologically</i> possible &#402;<sub>11</sub>: """
    print round(result["minhap"],5)
    print """<br>Maximum <i>biologically</i> possible &#402;<sub>11</sub>: """
    print round(result["maxhap"],5)
        
    print """</p><p><b>Number of impossible solutions: </b>"""
    print str(3-(result["alphaposs"] + result["betaposs"] + result["gammaposs"]))
    if result["alphaposs"] == 0:
        print """<br><b>&#945;: </b>&#402;<sub>11</sub> = """
        print round(result["alpha"],5)
    if result["betaposs"] == 0:
        print """<br><b>&#946;: </b>&#402;<sub>11</sub> = """
        if result["beta"] == "Not a real root":
            print result["beta"]
        else: 
            print round(result["beta"],5)
    if result["gammaposs"] == 0:
        print """<br><b>&#947;: </b>&#402;<sub>11</sub> = """
        if result["gamma"] == "Not a real root":
            print result["gamma"]
        else:
            print round(result["gamma"],5)
        
    print """</p><p><b>SNP 1 allele 1 frequency = </b> """
    print round(result["p"],3)
    print """<br><b>SNP 2 allele 1 frequency = </b> """
    print round(result["q"],3)
    print """</p><H3 style="border-top: 1px dashed rgb(187, 187, 187); background: rgb(255, 255, 255) none repeat scroll 0% 50%; -moz-background-clip: initial; -moz-background-origin: initial; -moz-background-inline-policy: initial; text-align: center; position: relative; width: 800px; z-index: 1;"><a name="note">Notes</a></H3>"""
    print """<ul><li style="width:700px;">Here &#402;<sub>11</sub> refers to an estimated haplotype frequency for allele 1 at locus 1 """
    print """and allele 1 at locus 2 (likewise for other haplotypes 12,21 and 22). The character &#402; should have a "hat" to indicate that it is estimated """
    print """but html limitations prevent this.</li>"""
    print """<li style="width:700px;">&#402;<sub>11</sub> is based on direct solution of the cubic equation expressing the phase uncertain """
    print """double heterozygotes (middle square of the 3x3) and overall model in terms of estimated &#402;<sub>11</sub> """
    print """assuming: <ol><li>random mating</li><li>Hardy-Weinberg equilibrium at both loci</li></ol></li>"""
    print """<li style="width:700px;">Expectation-maximisation algorithms rely on an iteration rather than direct solution. The direct """
    print """solution will display both the most likely (which EM should ideally reach) and other possible solutions """
    print """where iterations may converge in error.</li>"""
    print """<li style="width:700px;"><i>&#935;</i><sup>2</sup> values under the 3x3 represent difference between observed and expected diplotype frequencies. The number of degrees of freedom is equal to the number of observations (diplotype counts) minus four estimated parameters which are either three haplotypes (the fourth can be inferred) and D, or one haplotype, two allele frequencies and D. If nine different diplotypes are observed the number of degrees of freedom is therefore five. For each empty cell in the 3x3 the number of degrees of freedom is reduced by one. If the user knows there are only three haplotypes present (and therefore six diplotypes) then there are only three estimated parameters (D is inferred by the three haplotype frequencies) and 3df. It is important to note that in the latter case neither cubic solution nor iteration is necessary as the haplotype frequencies can be directly counted from the diplotype data. If the user believes that there are only three alleles and hence six diplotypes, but there are non-zero values for any of the other three possible diplotypes, then reconsideration of the technical veracity of the data and of the homogeneity of the population sample would be wise.</li>"""
    print """<li style="width:700px;">In perfect Hardy-Weinberg equilibrium "expected" numbers will be whole numbers matching the "observed" diplotype numbers. However, imperfect Hardy-Weinberg proportions will result in impossible "fractions" of individuals, which are shown for comparison to observations."""
    print """</ul>"""
    print """<hr>"""
    print """<a href="/software/cubex/">Return to input</a>"""
    print """<hr>"""
print """</div></body></html>"""
referrer = "None"
ipaddress = "None"
remotehost = "None"
querystring = "None"
for k,v in os.environ.items():
    if str(k.strip()) == "REMOTE_ADDR":
        ipaddress = str(v)
        thistime = time.time()
        thistimelocal = time.ctime(thistime)
    if str(k.strip()) == "HTTP_REFERER":
        referrer = str(v)
    if str(k.strip()) == "REMOTE_HOST":
        remotehost = str(v)
    if str(k.strip()) == "REQUEST_URI":
        querystring = str(v)
fileHandle = open ("../data/cubexusestats.txt", "a")
fileHandle.write(str(querystring) + "\t" + str(thistimelocal) + "\t" + 
str(ipaddress) + "\t" + str(remotehost) + "\t" + str(referrer) + "\n")
fileHandle.close()
