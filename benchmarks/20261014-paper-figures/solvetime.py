# Arithmetic behind paper/rho-blackwell Section 16 (Table "solvetime"): expected work 2^60.9,
# rho spread exp(-pi c^2/4), cutoff-32 detection delay, Modal list cost from the
# survey's iterations per dollar, and the NIST P-256/P-384 expected work sqrt(pi n/4).
# Rates are the recorded ones quoted in the paper; nothing is measured here.
import math
W=2**60.9; YR=365.25*86400; DPI=2**28.41
print("W=%.4g  W/yr=%.3g B/s"%(W,W/YR/1e9))
q=lambda p: math.sqrt(-4*math.log(p)/math.pi)
print("median c=%.3f p90=%.3f p99=%.3f"%(q(.5),q(.1),q(.01)))
# conditional remaining given c0 done
c0=2**56.68/W
rem=math.exp(math.pi*c0*c0/4)*math.erfc(c0*math.sqrt(math.pi)/2)
print("c0=%.4f remaining expected=%.4f"%(c0,rem))
# survey rows: name, rate, T/$ (survey), survey rate
rows=[("RTX PRO 6000 campaign",14.64,17.95,15.116),("RTX PRO 6000 best sigma",15.884,17.95,15.116),
("RTX PRO 4500 bridge",12.499/1.006696,None,None),("RTX PRO 4500 sigma",5.107,None,None),
("B200",8.840,5.09,8.840),("B300",7.575,3.84,7.575),("H200",7.632,6.05,7.632),("H100",7.540,6.87,7.540),
("L40S",8.696,16.04,8.696),("L4",2.561,11.54,2.561),("A100-80",4.632,6.67,4.632),("A10",2.423,7.92,2.423),("T4",0.541,3.30,0.541),
("F2 FPGA",8.25,None,None),("ceiling 24",24,None,None),("ceiling 26",26,None,None),("GTX295 2010",0.064,None,None),("CPU core AVX512",0.0126,None,None)]
for n,r,tpd,sr in rows:
    yrs=W/(r*1e9)/YR; hrs=W/(r*1e9)/3600
    s="%-26s %7.3f  years %8.2f  hours %9.0f  d@100 %7.2f d@1000 %6.2f"%(n,r,yrs,hrs,hrs/100/24,hrs/1000/24)
    if tpd:
        price=sr*1e9*3600/(tpd*1e12); s+="  $/h %.2f cost $%.0fk"%(price,hrs*price/1e3)
    print(s)
# latency on 6000
for r in (14.64,15.884,14.11):
    lat=6160384*DPI/(r*1e9); print("latency at %.2f: %.1f h"%(r,lat/3600))
print("w34 latency", 6160384*2**25.27/14.64e9/3600)
# number of 6000s for 1 week/1 month/1 year, at 15.884, mean
for d in (7,30,365):
    print(d,"days ->", W/(15.884e9*d*86400),"GPUs (no latency)", W/(14.64e9*d*86400))
# NIST
for name,n,r,b in (("P256",0xFFFFFFFF00000000FFFFFFFFFFFFFFFFBCE6FAADA7179E84F3B9CAC2FC632551,1.632,12.2),
                ("P384",int("FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFC7634D81F4372DDF581A0DB248B0A77AECEC196ACCC52973",16),0.548,None)):
    it=math.sqrt(math.pi*n/4); s=it/(r*1e9)
    print(name,"iters %.3g = 2^%.2f  sec %.3g years %.3g"%(it,math.log2(it),s,s/YR), "bound yrs %.3g"%(it/(b*1e9)/YR) if b else "")
    print("  on 1e9 GPUs years %.3g; vs ecc2k130 ratio of work %.3g"%(s/YR/1e9, it/W))
