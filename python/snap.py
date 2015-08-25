import matplotlib.pyplot as plt
import numpy 		 as np

Npart = 1024
NN    = Npart+2
index = [ 16, 27, 29, 52, 85, 119, 121, 123, 138, 147, 149, 186, 267, 276, 288, 294, 317, 324, 344, 438, 441, 564, 565, 568, 596, 601,
          602,811,924,927,1241,1242, 1342, 1564, 1847, 1851, 1852, 1853, 1858, 1866, 1877, 1888]

counter = 0
c2    = 0
f = open("xyz.xyz")

flag = False
#conf = [][] ,527,532,534,537,538,539,542,547
confx = []
confy = []
x = []
y = []
for line in f:
    counter = counter + 1
    ###if counter < 3: print line
    ###if c2%Npart == 0: c2 = 0
    if counter%NN == 3:
       a = counter/NN
       if a in index: 
           flag = True
       else: flag = False
    if flag and (counter%NN == 0 or counter%NN > 2):
       temp = line.split()
       #print temp
       x.append(float(temp[1]))
       y.append(float(temp[2]))
       c2 = c2 + 1
       if c2%Npart == 0: 
           confx.append(x)
           confy.append(y)
           c2 = 0
           x  = []
           y  = []


frame = len(confx)
print len(confx)

c3 = 0

while (c3 < frame):
    f, (ax1, ax2, ax3) = plt.subplots(3, sharex=True, sharey=True)
    ax1.plot(confx[c3], confy[c3], 'g-', lw=2)
    ax1.xaxis.set_visible(False)
    ax1.yaxis.set_visible(False)
    ax1.set_aspect('equal')
    ax1.set_ylim([0, 36])

    ax2.plot(confx[c3+1], confy[c3+1], 'g-', lw=2)
    ax2.xaxis.set_visible(False)
    ax2.yaxis.set_visible(False)
    ax2.set_aspect('equal')
    ax2.set_ylim([0, 36])

    ax3.plot(confx[c3+2], confy[c3+2], 'g-', lw=2)
    #ax3.xaxis.set_visible(False)
    ax3.yaxis.set_visible(False)
    ax3.set_aspect('equal')
    ax3.set_ylim([0, 36])


    #ax1.set_frame_on(False)
    #ax2.set_frame_on(False)
    #ax3.set_frame_on(False)
    #ax1.set_title(str(index[c3]))
    #ax2.set_title(str(index[c3+1]))
    #ax3.set_title(str(index[c3+2]))

    ax1.annotate(str(index[c3]), xy=(.1, 0.87),  xycoords='axes fraction',
                horizontalalignment='center', verticalalignment='center')
    ax2.annotate(str(index[c3+1]), xy=(.1, 0.87),  xycoords='axes fraction',
                horizontalalignment='center', verticalalignment='center')
    ax3.annotate(str(index[c3+2]), xy=(.1, 0.87),  xycoords='axes fraction',
                horizontalalignment='center', verticalalignment='center')
    ax1.axhline(y=0.000, xmin=0, xmax=1, c='r', lw=3)
    ax1.axhline(y=36.000, xmin=0, xmax=1, c='r', lw=3)
    ax2.axhline(y=0.000, xmin=0, xmax=1, c='r', lw=3)
    ax2.axhline(y=36.000, xmin=0, xmax=1, c='r', lw=3)
    ax3.axhline(y=0.000, xmin=0, xmax=1, c='r', lw=3)
    ax3.axhline(y=36.000, xmin=0, xmax=1, c='r', lw=3)

    f.subplots_adjust(hspace=0.15)
#plt.axis('off')
#plt.gca().set_aspect('equal')
    #plt.show()
    f.savefig(str(c3/3)+".eps", format="eps")
    c3 = c3+3

#f, (ax1, ax2, ax3) = plt.subplots(3, sharex=True, sharey=True)
#ax1.plot(x, y, 'g-', lw=3)


'''        
plt.plot(x, y, 'g-', lw=3)
plt.axhline(y=0.000, xmin=0, xmax=1, c='r', lw=5)
plt.axhline(y=80.000, xmin=0, xmax=1, c='r', lw=5)

#plt.axes(frameon=False)
plt.axis('off')
plt.gca().set_aspect('equal')


#plt.ylim(3,9)
#plt.xlim(3,95)

#plt.ylabel('bubble length', fontsize=25)
#plt.xlabel('sequence', fontsize=25)
#plt.tick_params( labelsize=20)
#plt.title('$\\lambda$ phage sequence (33.0k-33.1k) at 400K', fontsize=26)

plt.show()
'''
    
