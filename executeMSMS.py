from Peptide import *
import sys
import gzip
import xml.etree.ElementTree as ET
from base64 import b64decode
from array import array
import matplotlib.pyplot as plt
if len(sys.argv)<4:
	print "You seem to be missing either your file name, your peptide sequence, or your scan number? Please submit all three in that order."
	sys.exit(1)


inputseq=sys.argv[2]
try:
	if "." in inputseq:
		print "Invalid peptide input. Filename found where peptide sequence should be."
		sys.exit(1)
	int(inputseq)
	print "Invalid peptide input. Scan number found where peptide sequence should be."
	sys.exit()
except ValueError:
	protein= Peptide(inputseq)
try:
	scannum=str(int(sys.argv[3]))
except ValueError:
	print "Invalid Scan number input. Please input file name, peptide sequence, and scan number, in that order."
	sys.exit(1)
filename= sys.argv[1]
try:
	filenm= gzip.open(filename)
except IOError:
	print "Invalid Filename input; none found in directory. Please input file name, peptide sequence, and scan number, in that order."
	sys.exit(1)

#bion and yion calculation for peptide input
bion_dict={}
yion_dict={}
for i in range(1,protein.seqlen()+1):
	b_fragment= Peptide(protein.seq[0:i])
	bion_dict[b_fragment.seq]=b_fragment.b_ion()
for i in range(protein.seqlen()-1,-1,-1):
	y_fragment= Peptide(protein.seq[i:])
	yion_dict[y_fragment.seq]=y_fragment.y_ion()

# print "Original Protein sequence", protein.seq
# bionlist=sorted(bion_dict.items(), key = lambda p: p[1])
# print "sorted b_ionlist:\n", bionlist,"\n\n"
# yionlist=sorted(yion_dict.items(), key = lambda p: p[1])
# print "sorted y_ionlist:\n", yionlist,"\n\n"

						#******************#



ns= "{http://sashimi.sourceforge.net/schema/}"
for event,ele in ET.iterparse(filenm):
	if (ele.tag==ns+"scan") & ("num" in ele.attrib):
		if (ele.attrib['num']==scannum): 
			#print "ele.tag:", ele.tag
			peaklist= ele.find(ns+"peaks")
			# peaks elt is the XML element corresponding to the peaks list
			peaks = array('f',b64decode(peaklist.text))
			if sys.byteorder != 'big':
		  		peaks.byteswap()
				mzs = peaks[::2]
				ints = peaks[1::2]		
#mzs= [round(value,1) for value in mzs]
mzdict= dict(zip([round(value,1) for value in mzs],ints))

#create plot points and labels
bionlength_labels=[]
yionlength_labels=[]
plotmzbion=[]
plotmzyion=[]
for mz in mzdict:
	for bion in bion_dict:
		if mz==bion_dict[bion]:
			bionlength_labels.append(len(bion))
			plotmzbion.append((mz,mzdict[mz]))
	for yion in yion_dict:
		if mz==yion_dict[yion]:
			yionlength_labels.append(len(yion))
			plotmzyion.append((mz,mzdict[mz]))
								#******************#
#plot data points.


fig=plt.figure()
ax=fig.add_subplot(111)
plt.stem(*zip(*plotmzbion),linefmt='b-',markerfmt="b ",label="b_ion",basefmt='b-')
for label, xy in zip(bionlength_labels, plotmzbion):
	ax.annotate('(b%s)' %label, xy=xy, textcoords='offset points',xytext=(-10,5),fontsize=10)
plt.stem(*zip(*plotmzyion),linefmt='r-',markerfmt="r ",label="y_ion",basefmt='r-')
for label, xy in zip(yionlength_labels, plotmzyion):
	ax.annotate('(y%s)' %label, xy=xy, textcoords='offset points',xytext=(-10,5),fontsize=10)
#set axes
usedints=[mzbion[1] for mzbion in plotmzbion]+[mzyion[1] for mzyion in plotmzyion]
usedmzs=[mzbion[0] for mzbion in plotmzbion]+[mzyion[0] for mzyion in plotmzyion]
plt.ylim([0,max(usedints)*1.15])
plt.xlim([min(usedmzs)*.85,max(usedmzs)*1.15])
plt.xlabel('M/Z')
plt.ylabel('Intensity')

#create annotated protein with b and y ion divisions
bseqlist= list(protein.seq)
yseqlist=list(protein.seq)
for i in range(protein.seqlen()):
	if i+1 in bionlength_labels:
		bseqlist[i]=bseqlist[i]+"}"
	if protein.seqlen()-i in yionlength_labels:
		yseqlist[i]="{"+yseqlist[i]
bseqlist= ''.join(bseqlist)
yseqlist= ''.join(yseqlist)
#show annotated protein sequence on figure
plt.legend()
ax.set_title(protein.seq+", Scan # "+sys.argv[3])
plt.text(0.008,0.970,bseqlist, verticalalignment='bottom',horizontalalignment='left',color='blue',transform=ax.transAxes)
plt.text(0.008,0.940,yseqlist, verticalalignment='bottom',horizontalalignment='left',color='red',transform=ax.transAxes)
plt.show()

