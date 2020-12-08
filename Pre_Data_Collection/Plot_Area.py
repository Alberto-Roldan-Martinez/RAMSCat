



import sys
import numpy as np
import matplotlib.pyplot as plt


#DataSet = np.loadtxt(sys.argv[1], comments="#")                    # import and reads data
#DataSet_2 = np.loadtxt(sys.argv[2], comments="#")
DataFile = open(sys.argv[1])
DataSet = DataFile.readlines()[1:]
DataSet = [DataSet[i].split() for i in range(len(DataSet)-1) if float(DataSet[i][0]) > 0]


def annotation(note,x0,y0,x,y):
    return  plt.annotate(note, xy=(x0,y0), xycoords="data", size=14,
                               xytext=(x,y), textcoords="data",
                               arrowprops=dict(arrowstyle="<-", color="b", fc="0.6"),
                               horizontalalignment="right", verticalalignment="top")

def notes(i):
    line = DataSet[i]
    if line[4] != str(0):
        note = str("$\it p$("+line[2]+")$\cdot$("+line[3]+")"+line[4])
    else:
        note = str("$\it p$("+line[2]+")$\cdot$("+line[3]+")")
    return note


# What to plot?
X = [float(DataSet[i][0]) for i in range(len(DataSet)) if float(DataSet[i][0]) > 0]
Y = [float(DataSet[i][1]) for i in range(len(DataSet)) if float(DataSet[i][0]) > 0]

plt.plot(np.linspace(0,max(X)+10,2),np.linspace(0,max(X)+10,2), "k-", lw=1.5)
plt.plot(X, Y,"rs", markersize=7)
#plt.plot(X2, Y2,"rs", fillstyle = "none", markersize=7)


plt.xlabel("Area /$\\AA ^{2}$", fontsize=16)
plt.ylabel("Predicted Area /$\\AA ^{2}$", fontsize=16)
plt.tick_params(axis='both',labelrotation=0,labelsize=14)           # custimise tick labels

i = 0  ;  x = X[i]+50  ;  y = Y[i]     ; annotation(notes(i),X[i],Y[i],x,y)
i = 1  ;  x = X[i]+50  ;  y = Y[i]+5   ; annotation(notes(i),X[i],Y[i],x,y)
i = 2  ;  x = X[i]+75  ;  y = Y[i]     ; annotation(notes(i),X[i],Y[i],x,y)
i = 3  ;  x = X[i]+75  ;  y = Y[i]     ; annotation(notes(i),X[i],Y[i],x,y)
i = 4  ;  x = X[i]+85  ;  y = Y[i]     ; annotation(notes(i),X[i],Y[i],x,y)
i = 5  ;  x = X[i]+60  ;  y = Y[i]+5   ; annotation(notes(i),X[i],Y[i],x,y)
i = 6  ;  x = X[i]-20  ;  y = Y[i]+10  ; annotation(notes(i),X[i],Y[i],x,y)
i = 7  ;  x = X[i]-15  ;  y = Y[i]+15  ; annotation(notes(i),X[i],Y[i],x,y)
i = 8  ;  x = X[i]+70  ;  y = Y[i]     ; annotation(notes(i),X[i],Y[i],x,y)
i = 9  ;  x = X[i]-15  ;  y = Y[i]+3   ; annotation(notes(i),X[i],Y[i],x,y)
i = 10 ;  x = X[i]+140 ;  y = Y[i]+10  ; annotation(notes(i),X[i],Y[i],x,y)
i = 11 ;  x = X[i]-35  ;  y = Y[i]+20  ; annotation(notes(i),X[i],Y[i],x,y)
i = 12 ;  x = X[i]+60  ;  y = Y[i]+5   ; annotation(notes(i),X[i],Y[i],x,y)
i = 13 ;  x = X[i]+80  ;  y = Y[i]     ; annotation(notes(i),X[i],Y[i],x,y)
i = 14 ;  x = X[i]+120 ;  y = Y[i]+3   ; annotation(notes(i),X[i],Y[i],x,y)
i = 15 ;  x = X[i]+60  ;  y = Y[i]+8   ; annotation(notes(i),X[i],Y[i],x,y)
i = 16 ;  x = X[i]+100 ;  y = Y[i]+3   ; annotation(notes(i),X[i],Y[i],x,y)
i = 17 ;  x = X[i]+50  ;  y = Y[i]+3   ; annotation(notes(i),X[i],Y[i],x,y)
i = 18 ;  x = X[i]-5   ;  y = Y[i]+3   ; annotation(notes(i),X[i],Y[i],x,y)
i = 19 ;  x = X[i]-30  ;  y = Y[i]+20  ; annotation(notes(i),X[i],Y[i],x,y)
i = 20 ;  x = X[i]-10  ;  y = Y[i]+2   ; annotation(notes(i),X[i],Y[i],x,y)
i = 21 ;  x = X[i]-30  ;  y = Y[i]+5   ; annotation(notes(i),X[i],Y[i],x,y)
i = 22 ;  x = X[i]+20  ;  y = Y[i]+30  ; annotation(notes(i),X[i],Y[i],x,y)
i = 23 ;  x = X[i]-25  ;  y = Y[i]+30  ; annotation(notes(i),X[i],Y[i],x,y)
i = 24 ;  x = X[i]-20  ;  y = Y[i]+33  ; annotation(notes(i),X[i],Y[i],x,y)
i = 25 ;  x = X[i]-15  ;  y = Y[i]+40  ; annotation(notes(i),X[i],Y[i],x,y)
i = 26 ;  x = X[i]-18  ;  y = Y[i]+23  ; annotation(notes(i),X[i],Y[i],x,y)
i = 27 ;  x = X[i]+80  ;  y = Y[i]+3   ; annotation(notes(i),X[i],Y[i],x,y)
i = 28 ;  x = X[i]+50  ;  y = Y[i]+5   ; annotation(notes(i),X[i],Y[i],x,y)
i = 29 ;  x = X[i]-20  ;  y = Y[i]+22  ; annotation(notes(i),X[i],Y[i],x,y)
i = 30 ;  x = X[i]-20  ;  y = Y[i]     ; annotation(notes(i),X[i],Y[i],x,y)
i = 31 ;  x = X[i]+55  ;  y = Y[i]+3   ; annotation(notes(i),X[i],Y[i],x,y)
i = 32 ;  x = X[i]-20  ;  y = Y[i]+20  ; annotation(notes(i),X[i],Y[i],x,y)


plt.tight_layout()
plt.grid(True)
plt.savefig("Area.png", figsize=(11.69,16.53), clear=True, dpi=300, orientation='portrait',transparent=True)
plt.show()

