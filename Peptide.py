import sys

class Seq:
	def __init__(self,seq):
		self.seq=seq
	def seqlen(self):
		return len(self.seq)
class Peptide(Seq):
	mw={'A':71.04 ,'C':103.01, 'D':115.03, 'E':129.04, 'F':147.07,'G':57.02,'H':137.06,'I':113.08,'K':128.09,'L':113.08,'M':131.04,'N':114.04, 'P':97.05,'Q':128.06,'R':156.10,'S':87.03,'T':101.05,'V':99.07,'W':186.08,'Y':163.06}
	#valid_sym='ACDEFGHIKLMNPQRSTVWY'
	def molWt(self):
		try:
			return sum(map(self.mw.get,self.seq))
		except TypeError:
			print "Invalid Peptide sequence found while calculating b and y ions.",self.seq
			sys.exit(1)
	def b_ion(self):
		return round(self.molWt()+1,1)
	def y_ion(self):
		return round(self.molWt()+19,1)